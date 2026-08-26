#!/usr/bin/env Rscript

# Normalize the configured cohort-level testing-unit annotations.

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
        description = "Normalize cohort-level testing-unit annotations."
    )
    parser$add_argument("--input", required = TRUE)
    parser$add_argument("--output", required = TRUE)
    parser$add_argument("--threads", type = "integer", default = 1L)

    parsed <- parser$parse_args(args = args)
    require_positive_threads(parsed$threads)
    parsed
}


clean_unit_value <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x == "" | x == "." | x == "-"] <- NA_character_
    x
}


expand_units <- function(annotations) {
    if (nrow(annotations) == 0L) {
        annotations[, unit := character()]
        return(annotations)
    }

    pieces <- strsplit(
        clean_unit_value(annotations$unit),
        split = ",",
        fixed = TRUE
    )
    row_index <- rep.int(seq_len(nrow(annotations)), lengths(pieces))
    expanded <- annotations[row_index]
    expanded[, unit := clean_unit_value(unlist(pieces, use.names = FALSE))]
    expanded[!is.na(unit)]
}


prepare_testing_units <- function(annotations) {
    require_columns(
        annotations,
        c("CHROM", "POS", "REF", "ALT", "CADD_PHRED", "unit"),
        "Configured testing-unit annotation table"
    )

    annotations <- add_vcf_variant_identity(annotations)
    annotations <- expand_units(annotations)
    annotations <- unique(annotations, by = c("Key", "unit"))
    annotations[, .(Key, CADD_PHRED, unit)]
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    message("Reading configured testing-unit annotations: ", config$input)
    annotations <- fread(
        config$input,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )
    input_rows <- nrow(annotations)
    output <- prepare_testing_units(annotations)

    fwrite(
        output,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Input annotation rows: ", input_rows)
    message("Testing-unit rows written: ", nrow(output))
    message("Unique variants written: ", uniqueN(output$Key))
    message("Unique testing units written: ", uniqueN(output$unit))
    message("Wrote testing-unit table: ", config$output)
}


main()
