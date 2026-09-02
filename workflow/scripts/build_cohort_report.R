#!/usr/bin/env Rscript

# Snakemake adapter for cohort-level candidate/control reporting.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) source(file.path(snakemake@scriptdir, "bootstrap.R"))
options(scipen = 999)


workflow_config <- function() {
  list(
    candidates = as.character(snakemake@input[["candidates"]]),
    controls = snakemake@input[["controls"]],
    output_candidates = snakemake@output[["candidates"]],
    output_summary = snakemake@output[["summary"]],
    threads = snakemake@threads
  )
}


if (exists("snakemake", inherits = FALSE)) run_cohort_reporting(workflow_config())
