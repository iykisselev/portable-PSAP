#!/usr/bin/env Rscript

# Snakemake adapter for affected-sample candidate construction.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) source(file.path(snakemake@scriptdir, "bootstrap.R"))
options(scipen = 999)


workflow_config <- function() {
  pid <- as.character(snakemake@params[["pid"]])
  if (!nzchar(trimws(pid))) abort("Sample ID must be non-empty")
  list(
    input = snakemake@input[["psap"]], pid = pid,
    sex_value = parse_sex(as.character(snakemake@params[["sex"]])),
    candidate_id_separator = as.character(snakemake@params[["candidate_id_separator"]]),
    candidate_id_ordinal_width = as.integer(snakemake@params[["candidate_id_ordinal_width"]]),
    output = snakemake@output[["candidates"]], threads = snakemake@threads
  )
}


if (exists("snakemake", inherits = FALSE)) run_candidate_reporting(workflow_config())
