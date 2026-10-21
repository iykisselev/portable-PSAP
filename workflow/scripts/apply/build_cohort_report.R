#!/usr/bin/env Rscript

# Snakemake adapter for cohort-level candidate/control reporting.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}
options(scipen = 999)


workflow_config <- function() {
  list(
    candidates = as.character(snakemake@input[["candidates"]]),
    controls = snakemake@input[["controls"]],
    output_candidates = snakemake@output[["candidates"]],
    output_summary = snakemake@output[["summary"]],
    output_partition_summaries = as.character(snakemake@output[["partition_summaries"]]),
    partition_ids = as.character(snakemake@params[["partition_ids"]]),
    threads = snakemake@threads
  )
}


if (exists("snakemake", inherits = FALSE)) {
  config <- workflow_config()
  candidates <- data.table::rbindlist(lapply(config$candidates, read_tsv), use.names = TRUE, fill = TRUE)
  reports <- partitioned_cohort_reports(candidates, config$partition_ids)
  write_tsv(candidates, config$output_candidates, threads = config$threads)
  write_tsv(reports$main, config$output_summary, threads = config$threads)
  for (i in seq_along(config$partition_ids)) write_tsv(reports$partitions[[i]], config$output_partition_summaries[[i]], threads = config$threads)
}
