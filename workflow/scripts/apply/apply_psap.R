#!/usr/bin/env Rscript

# Apply precomputed PSAP lookup tables to one sample's annotated variants.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))
rm(script_path)

options(scipen = 999)

SCORING_OUTPUT_COLUMNS <- c("Key", "unit", "Dz.Model", "PSAP")


empty_scoring_output <- function() {
    data.table(
        Key = character(),
        unit = character(),
        `Dz.Model` = character(),
        PSAP = numeric()
    )
}


parse_cli <- function(args) {
    parser <- ArgumentParser(
        description = paste(
            "Calculate per-sample PSAP candidates from lookup tables."
        )
    )

    parser$add_argument(
        "--input",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed PSAP-ready sample variant TSV."
    )

    parser$add_argument(
        "--lookup-het",
        required = TRUE,
        metavar = "FILE",
        help = "Heterozygous-model PSAP lookup table."
    )

    parser$add_argument(
        "--lookup-hom",
        required = TRUE,
        metavar = "FILE",
        help = "Homozygous-model PSAP lookup table."
    )

    parser$add_argument("--lookup-hem-chrx", default = "", metavar = "FILE")

    parser$add_argument(
        "--score-max",
        required = TRUE,
        type = "double",
        metavar = "N",
        help = "Maximum CADD score represented by the lookup grid."
    )

    parser$add_argument(
        "--hemizygote",
        required = TRUE,
        metavar = "BOOL",
        help = "Whether to apply the X-linked hemizygote model."
    )

    parser$add_argument(
        "--output",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed PSAP candidate output TSV."
    )

    parser$add_argument(
        "--threads",
        type = "integer",
        default = 1L,
        metavar = "N",
        help = "Number of data.table threads to use [default: %(default)s]."
    )

    parsed <- parser$parse_args(args = args)

    if (is.na(parsed$score_max) || parsed$score_max <= 0) {
        abort("--score-max must be positive")
    }

    require_positive_threads(parsed$threads)

    parsed$hemizygote_model <- parse_bool(
        parsed$hemizygote,
        "--hemizygote"
    )

    if (parsed$hemizygote_model && !nzchar(trimws(parsed$lookup_hem_chrx))) {
        abort("--lookup-hem-chrx is required when the hemizygous model is enabled")
    }

    parsed$lookup_paths <- list(
        het = parsed$lookup_het,
        hom = parsed$lookup_hom,
        hem_chrX = parsed$lookup_hem_chrx
    )

    parsed
}


read_variants <- function(path) {
    variants <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    require_columns(
        variants,
        c("Key", "CHROM", "CADD_PHRED", "gt_class", "unit"),
        "PSAP-ready variant table"
    )

    variants[, CADD_PHRED := suppressWarnings(as.numeric(CADD_PHRED))]
    variants
}


read_lookup <- function(path) {
    lookup <- fread(
        path,
        sep = "\t",
        header = FALSE,
        quote = "",
        colClasses = "character",
        data.table = TRUE,
        showProgress = FALSE
    )

    if (ncol(lookup) < 2 || nrow(lookup) == 0) {
        abort("PSAP lookup table must contain an ID column and score columns: %s", path)
    }

    ids <- trimws(as.character(lookup[[1]]))
    values <- lookup[, -1, with = FALSE]
    value_matrix <- matrix(
        suppressWarnings(as.numeric(unlist(values, use.names = FALSE))),
        nrow = nrow(values),
        ncol = ncol(values)
    )

    # Accept an optional header row while keeping the normal headerless format.
    invalid_rows <- apply(value_matrix, 1, function(row) any(is.na(row)))
    if (invalid_rows[1] && any(!invalid_rows[-1])) {
        ids <- ids[-1]
        value_matrix <- value_matrix[-1, , drop = FALSE]
    }

    if (anyNA(value_matrix)) {
        abort("PSAP lookup table contains non-numeric score values: %s", path)
    }

    if (anyDuplicated(ids)) {
        abort("PSAP lookup table contains duplicate testing-unit IDs: %s", path)
    }

    list(
        ids = ids,
        values = value_matrix
    )
}


lookup_psap <- function(variants, lookup, score_max, score_column = "CADD_PHRED") {
    score <- pmin(
        pmax(variants[[score_column]], 0),
        score_max
    )
    score_grid <- seq(
        0,
        score_max,
        length.out = ncol(lookup$values)
    )
    score_bin <- findInterval(score, score_grid)
    score_bin[score_bin < 1] <- 1
    score_bin[score_bin > ncol(lookup$values)] <- ncol(lookup$values)

    row_index <- match(variants$.testing_unit_id, lookup$ids)

    if (anyNA(row_index)) {
        missing_ids <- head(
            unique(variants$.testing_unit_id[is.na(row_index)]),
            5
        )
        abort(
            "PSAP lookup table is missing testing-unit ID(s): %s",
            paste(missing_ids, collapse = ", ")
        )
    }

    lookup$values[cbind(row_index, score_bin)]
}


select_highest <- function(variants) {
    setorderv(
        variants,
        c(".testing_unit_id", "CADD_PHRED", "Key"),
        c(1, -1, 1),
        na.last = TRUE
    )

    variants[, head(.SD, 1), by = .testing_unit_id]
}


score_model <- function(
    variants,
    model,
    genotype_class,
    lookup,
    score_max
) {
    selected <- variants[
        gt_class == genotype_class &
            !is.na(CADD_PHRED) &
            !is.na(.testing_unit_id) &
            nzchar(.testing_unit_id)
    ]

    if (nrow(selected) == 0) {
        return(NULL)
    }

    require_lookup_units(selected, lookup, model)

    selected <- select_highest(selected)
    selected[, PSAP := lookup_psap(selected, lookup, score_max)]

    selected[
        ,
        .(
            Key = as.character(Key),
            unit = as.character(.testing_unit_id),
            `Dz.Model` = as.character(model),
            PSAP = as.numeric(PSAP)
        )
    ]
}


score_variants <- function(
    variants,
    lookup_paths,
    score_max,
    hemizygote_model
) {
    variants <- copy(variants)
    variants[, .testing_unit_id := clean_identifier(unit)]
    variants[, gt_class := tolower(trimws(as.character(gt_class)))]

    definitions <- DISEASE_MODELS[DISEASE_MODELS$lookup != "hem_chrX", ]

    lookup_names <- unique(definitions$lookup)
    lookups <- setNames(
        lapply(
            lookup_names,
            function(name) read_lookup(lookup_paths[[name]])
        ),
        lookup_names
    )

    results <- lapply(seq_len(nrow(definitions)), function(index) {
        definition <- definitions[index, ]
        score_model(
            variants,
            definition$model,
            definition$genotype_class,
            lookups[[definition$lookup]],
            score_max
        )
    })

    if (hemizygote_model) {
        hem_variants <- copy(variants)
        hem_variants[
            ,
            .chromosome := sub(
                "^CHR",
                "",
                toupper(trimws(as.character(CHROM)))
            )
        ]
        hem_definition <- DISEASE_MODELS[
            DISEASE_MODELS$lookup == "hem_chrX",
            ,
            drop = FALSE
        ]
        hem_lookup <- read_lookup(lookup_paths[[hem_definition$lookup]])

        results <- c(
            results,
            list(score_model(
                hem_variants[
                    .chromosome == "X" & gt_class == "hem"
                ],
                hem_definition$model,
                hem_definition$genotype_class,
                hem_lookup,
                score_max
            ))
        )
    }

    results <- Filter(Negate(is.null), results)

    if (length(results) == 0) {
        return(empty_scoring_output())
    }

    output <- rbindlist(results, use.names = TRUE, fill = FALSE)
    setcolorder(output, SCORING_OUTPUT_COLUMNS)
    output
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    message("Reading PSAP-ready variants: ", config$input)
    variants <- read_variants(config$input)
    message("Input rows: ", nrow(variants))

    output <- score_variants(
        variants,
        lookup_paths = config$lookup_paths,
        score_max = config$score_max,
        hemizygote_model = config$hemizygote_model
    )

    if (nrow(output) > 0) {
        setorderv(output, "PSAP", 1L, na.last = TRUE)
    }

    if (!identical(names(output), SCORING_OUTPUT_COLUMNS)) {
        abort(
            "Internal scoring output schema mismatch: expected %s; found %s",
            paste(SCORING_OUTPUT_COLUMNS, collapse = ", "),
            paste(names(output), collapse = ", ")
        )
    }

    fwrite(
        output,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Candidate rows written: ", nrow(output))
    message("Wrote PSAP results: ", config$output)
}


if (sys.nframe() == 0) {
    main()
}
