#!/usr/bin/env Rscript
# =======================================================================
# Apply sex-aware genotype classification to standardized genotype calls.
# =======================================================================

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
            "Apply sex-aware classification to standardized genotype calls."
        )
    )

    parser$add_argument(
        "--input",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed standardized genotype TSV."
    )

    parser$add_argument(
        "--output",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed classified genotype TSV output."
    )

    parser$add_argument(
        "--sex",
        required = TRUE,
        metavar = "N",
        help = "Numeric sex code: 1 male, 2 female, other or missing unknown."
    )

    parser$add_argument(
        "--hemizygote-model",
        required = TRUE,
        metavar = "BOOL",
        help = "Whether to classify male X/Y calls as hemizygous: True or False."
    )

    parser$add_argument("--threads", type = "integer", default = 1L)

    parsed <- parser$parse_args(args = args)
    require_positive_threads(parsed$threads)
    parsed$sex_value <- parse_sex(parsed$sex)
    parsed$hemizygote <- parse_bool(
        parsed$hemizygote_model,
        "--hemizygote-model"
    )

    parsed
}


read_genotypes <- function(path) {
    genotypes <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    require_columns(
        genotypes,
        c("Key", "CHROM", "GT", "gt_class"),
        "Standardized genotype table"
    )

    genotypes
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    message("Reading standardized genotypes: ", config$input)
    genotypes <- read_genotypes(config$input)

    classified <- apply_sex_classification(
        genotypes,
        sex = config$sex_value,
        hemizygote_model = config$hemizygote
    )

    message(
        "Male X/Y heterozygous calls removed: ",
        classified$removed_male_xy_het
    )
    message("Output rows: ", nrow(classified$genotypes))

    fwrite(
        classified$genotypes,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Wrote classified genotypes: ", config$output)
}


main()
