#!/usr/bin/env Rscript

# =============================================================================
# Join one sample's standardized genotypes to one testing-unit annotation table.
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
            "Join standardized sample genotypes to testing-unit annotations."
        )
    )

    parser$add_argument(
        "--gt",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed standardized genotype TSV."
    )

    parser$add_argument(
        "--unit-table",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed testing-unit annotation TSV."
    )

    parser$add_argument(
        "--output",
        required = TRUE,
        metavar = "FILE",
        help = "Compressed joined TSV output path."
    )

    parser$add_argument(
        "--threads",
        type = "integer",
        default = 1,
        metavar = "N",
        help = "Number of data.table threads to use [default: %(default)s]."
    )

    parsed <- parser$parse_args(args = args)

    require_positive_threads(parsed$threads)

    parsed
}


read_table <- function(path, description) {
    table <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    if (ncol(table) == 0) {
        abort("%s is empty: %s", description, path)
    }

    table
}


join_annotations <- function(genotypes, unit_table) {
    annotation_columns <- c(
        "Key",
        "CADD_PHRED",
        "unit"
    )

    genotype_columns <- c(
        "Key",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "GT",
        "gt_class"
    )

    require_columns(
        unit_table,
        annotation_columns,
        "Testing-unit annotation table"
    )
    require_columns(
        genotypes,
        genotype_columns,
        "Standardized genotype table"
    )

    if (anyDuplicated(genotypes$Key)) {
        abort("Standardized genotype table contains duplicate Key values")
    }

    joined <- merge(
        unit_table[, annotation_columns, with = FALSE],
        genotypes[, genotype_columns, with = FALSE],
        by = "Key",
        all = FALSE,
        sort = FALSE
    )

    joined[, c(annotation_columns, setdiff(genotype_columns, "Key")), with = FALSE]
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    message("Reading testing-unit annotations: ", config$unit_table)
    unit_table <- read_table(
        config$unit_table,
        "Testing-unit annotation table"
    )

    message("Reading standardized genotypes: ", config$gt)
    genotypes <- read_table(
        config$gt,
        "Standardized genotype table"
    )

    joined <- join_annotations(
        genotypes,
        unit_table
    )

    fwrite(
        joined,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Testing unit: ", config$unit)
    message("Joined rows written: ", nrow(joined))
    message("Wrote joined annotations: ", config$output)
}


main()
