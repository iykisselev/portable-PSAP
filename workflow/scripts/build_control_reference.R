#!/usr/bin/env Rscript

# Build a control reference from every carried variant in unaffected samples.

suppressPackageStartupMessages({
    library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
    source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


workflow_config <- function() {
    inputs <- as.character(snakemake@input[["annotated"]])
    psap_inputs <- as.character(snakemake@input[["psap"]])
    pids <- as.character(snakemake@params[["pids"]])
    sexes <- as.character(snakemake@params[["sexes"]])
    lengths <- c(
        input = length(inputs), psap_input = length(psap_inputs),
        pid = length(pids), sex = length(sexes)
    )
    if (length(unique(lengths)) != 1L) {
        abort("Control input counts differ: %s", paste(names(lengths), lengths, sep = "=", collapse = ", "))
    }
    list(
        input = inputs,
        psap_input = psap_inputs,
        pid = pids,
        sex_value = vapply(sexes, parse_sex, numeric(1)),
        hemizygote_model = parse_bool(
            as.character(snakemake@params[["hemizygote"]]), "model.hemizygote"
        ),
        deduplicate = isTRUE(snakemake@params[["deduplicate"]]),
        output = snakemake@output[["control_variants"]],
        threads = snakemake@threads
    )
}


empty_reference <- function() {
    data.table(
        pid = character(),
        phenotype = character(),
        sex = numeric(),
        unit = character(),
        model = character(),
        variant_key = character(),
        PSAP = numeric()
    )
}


make_evidence <- function(variants, row_filter, model, unit = variants$unit) {
    selected <- variants[row_filter]
    if (nrow(selected) == 0L) return(NULL)

    selected[
        ,
        .(
            unit = unit[row_filter],
            model = model,
            variant_key
        )
    ]
}


derive_control_evidence <- function(
    path,
    pid,
    sex,
    hemizygote_model
) {
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
        c("Key", "gt_class"),
        sprintf("Control annotation table for %s", pid)
    )

    variants[, variant_key := clean_identifier(Key)]
    variants[, gt_class := tolower(clean_identifier(gt_class))]
    require_columns(
        variants,
        "unit",
        sprintf("Control annotation table for %s", pid)
    )
    variants[, unit := clean_identifier(unit)]

    invalid <- is.na(variants$variant_key) | is.na(variants$gt_class) | is.na(variants$unit)
    if (any(invalid)) {
        abort(
            "Control %s contains %d row(s) with missing key, genotype class, or unit",
            pid,
            sum(invalid)
        )
    }

    definitions <- DISEASE_MODELS[DISEASE_MODELS$lookup != "hem_chrX", ]

    evidence <- lapply(seq_len(nrow(definitions)), function(index) {
        definition <- definitions[index, ]
        make_evidence(
            variants,
            variants$gt_class == definition$genotype_class,
            definition$model
        )
    })

    if (hemizygote_model && any(variants$gt_class == "hem")) {
        require_columns(
            variants,
            "CHROM",
            sprintf("Hemizygous control annotation table for %s", pid)
        )
        chromosome <- sub(
            "^CHR",
            "",
            toupper(trimws(as.character(variants$CHROM)))
        )
        hem_definition <- DISEASE_MODELS[
            DISEASE_MODELS$lookup == "hem_chrX",
            ,
            drop = FALSE
        ]
        evidence <- c(
            evidence,
            list(make_evidence(
                variants,
                variants$gt_class == hem_definition$genotype_class &
                    chromosome == "X",
                hem_definition$model
            ))
        )
    }

    evidence <- Filter(Negate(is.null), evidence)
    if (length(evidence) == 0L) return(empty_reference())

    output <- unique(
        rbindlist(evidence, use.names = TRUE),
        by = c("unit", "model", "variant_key")
    )
    output[
        ,
        `:=`(
            pid = as.character(pid),
            phenotype = "unaffected",
            sex = sex,
            PSAP = NA_real_
        )
    ]
    setcolorder(
        output,
        c("pid", "phenotype", "sex", "unit", "model", "variant_key", "PSAP")
    )
    output
}


read_psap_scores <- function(path, pid) {
    scores <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    require_columns(
        scores,
        c("Key", "unit", "Dz.Model", "PSAP"),
        sprintf("Control PSAP table for %s", pid)
    )

    scores[, unit := clean_identifier(unit)]
    scores[, model := clean_identifier(`Dz.Model`)]
    scores[, variant_key := clean_identifier(Key)]
    scores[, PSAP := suppressWarnings(as.numeric(PSAP))]
    scores <- scores[
        !is.na(unit) & !is.na(model) & !is.na(variant_key),
        .(unit, model, variant_key, PSAP)
    ]

    scores[
        ,
        .(PSAP = if (all(is.na(PSAP))) NA_real_ else min(PSAP, na.rm = TRUE)),
        by = .(unit, model, variant_key)
    ]
}


attach_psap <- function(evidence, scores) {
    if (nrow(evidence) == 0L || nrow(scores) == 0L) return(evidence)

    evidence[, PSAP := NULL]
    scores[evidence, on = c("unit", "model", "variant_key")][
        ,
        .(pid, phenotype, sex, unit, model, variant_key, PSAP)
    ]
}


run <- function(config) {
    setDTthreads(config$threads)

    if (length(config$input) == 0L) {
        reference <- empty_reference()
    } else {
        controls <- vector("list", length(config$input))

        for (i in seq_along(config$input)) {
            evidence <- derive_control_evidence(
                config$input[[i]],
                config$pid[[i]],
                config$sex_value[[i]],
                config$hemizygote_model
            )
            scores <- read_psap_scores(config$psap_input[[i]], config$pid[[i]])
            controls[[i]] <- attach_psap(evidence, scores)
        }

        reference <- rbindlist(controls, use.names = TRUE, fill = FALSE)
    }

    if (config$deduplicate && nrow(reference) > 0L) {
        before <- nrow(reference)
        reference <- unique(reference, by = names(reference))
        message("Exact duplicate rows removed: ", before - nrow(reference))
    }

    if (nrow(reference) > 0L) {
        setorderv(
            reference,
            c("unit", "model", "variant_key", "pid"),
            c(1L, 1L, 1L, 1L)
        )
    }

    fwrite(
        reference,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Unaffected individuals: ", length(config$input))
    message("Control evidence rows written: ", nrow(reference))
    message("Wrote control reference: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
    run(workflow_config())
}
