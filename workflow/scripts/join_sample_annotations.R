#!/usr/bin/env Rscript

# =============================================================================
# Join one sample's standardized genotypes to one testing-unit annotation table.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
    source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


workflow_config <- function() {
    list(
        gt = snakemake@input[["gt"]],
        unit_table = snakemake@input[["unit_tbl"]],
        output = snakemake@output[["annotated_gt"]],
        threads = snakemake@threads
    )
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

    joined <- merge(
        unit_table[, annotation_columns, with = FALSE],
        genotypes[, genotype_columns, with = FALSE],
        by = "Key",
        all = FALSE,
        sort = FALSE
    )

    joined[, c(annotation_columns, setdiff(genotype_columns, "Key")), with = FALSE]
}


run <- function(config) {
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

    message("Joined rows written: ", nrow(joined))
    message("Wrote joined annotations: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
    run(workflow_config())
}
