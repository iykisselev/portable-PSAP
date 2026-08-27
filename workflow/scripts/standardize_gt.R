#!/usr/bin/env Rscript

# =============================================================================
# Standardize one-sample genotype calls for downstream PSAP sample preparation.
# =============================================================================

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))
rm(script_path)

options(scipen = 999)


parse_cli <- function(args) {
    parser <- ArgumentParser(
        description = paste(
            "Standardize one-sample genotype calls from a bcftools query TSV."
        )
    )

    parser$add_argument(
        "--input",
        required = TRUE,
        metavar = "FILE",
        help = "Gzipped TSV containing CHROM, POS, REF, ALT, and GT fields."
    )

    parser$add_argument(
        "--output",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed TSV output path."
    )

    parser$add_argument(
        "--hemizygote-model",
        required = TRUE,
        metavar = "BOOL",
        help = "Whether haploid X/Y calls are supported: True or False."
    )

    parser$add_argument("--threads", type = "integer", default = 1L)

    parsed <- parser$parse_args(args = args)
    require_positive_threads(parsed$threads)
    parsed$hemizygote <- parse_bool(
        parsed$hemizygote_model,
        "--hemizygote-model"
    )
    parsed
}


read_genotypes <- function(path) {
    columns <- c("CHROM", "POS", "REF", "ALT", "GT")

    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    } else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)

    if (length(readLines(connection, n = 1L, warn = FALSE)) == 0L) {
        return(data.table(
            CHROM = character(),
            POS = integer(),
            REF = character(),
            ALT = character(),
            GT = character()
        ))
    }

    genotypes <- fread(
        path,
        sep = "\t",
        header = FALSE,
        col.names = columns,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    if (ncol(genotypes) != length(columns)) {
        abort(
            "Expected five input columns: %s; found %d",
            paste(columns, collapse = ", "),
            ncol(genotypes)
        )
    }

    genotypes
}


deduplicate_by_key <- function(genotypes) {
    if (nrow(genotypes) == 0L) {
        return(genotypes)
    }

    calls_by_key <- genotypes[
        ,
        .(
            n_calls = uniqueN(GT),
            calls = paste(sort(unique(GT)), collapse = ",")
        ),
        by = Key
    ]
    conflicts <- calls_by_key[n_calls > 1L]

    if (nrow(conflicts) > 0L) {
        examples <- head(
            paste0(conflicts$Key, " [", conflicts$calls, "]"),
            5L
        )
        abort(
            "Conflicting non-reference calls for %d variant key(s): %s",
            nrow(conflicts),
            paste(examples, collapse = "; ")
        )
    }

    unique(genotypes, by = "Key")
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    message("Reading genotype calls: ", config$input)
    genotypes <- read_genotypes(config$input)
    message("Input rows: ", nrow(genotypes))

    genotypes <- add_vcf_variant_identity(genotypes)
    before <- nrow(genotypes)
    genotypes <- standardize_genotypes(
        genotypes,
        hemizygote_model = config$hemizygote
    )
    message("Missing/reference genotype rows removed: ", before - nrow(genotypes))

    before <- nrow(genotypes)
    genotypes <- deduplicate_by_key(genotypes)
    message("Duplicate Key rows removed: ", before - nrow(genotypes))

    output_columns <- c(
        "Key",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "GT",
        "gt_class"
    )

    fwrite(
        genotypes[, output_columns, with = FALSE],
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Output rows: ", nrow(genotypes))
    message("Wrote standardized genotypes: ", config$output)
}


if (sys.nframe() == 0) {
    main()
}
