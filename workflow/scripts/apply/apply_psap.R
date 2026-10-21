#!/usr/bin/env Rscript

# Snakemake adapter for per-sample PSAP scoring.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}
options(scipen = 999)


workflow_config <- function() {
  hemizygote_model <- require_boolean(snakemake@params[["hemizygote"]], "model.hemizygote")
  score_max <- as.numeric(snakemake@params[["score_max"]])
  if (is.na(score_max) || score_max <= 0) abort("model.score_max must be positive")
  lookup_paths <- as.character(snakemake@input[["lookup_tables"]])
  lookup_keys <- as.character(snakemake@params[["lookup_keys"]])
  if (length(lookup_paths) != length(lookup_keys)) abort("Partition/model lookup paths do not match keys")
  names(lookup_paths) <- lookup_keys
  lookup_source <- resolve_lookup_contract(
    as.character(snakemake@params[["lookup_contract"]]),
    lookup_paths,
    expected_score_max = score_max,
    expected_bin_width = as.character(snakemake@params[["lookup_bin_width"]])
  )
  lookup_source$partition_ids <- as.character(snakemake@params[["partition_ids"]])
  list(input = snakemake@input[["annotated_gt"]], lookup_source = lookup_source,
       hemizygote_model = hemizygote_model, output = snakemake@output[["psap"]], threads = snakemake@threads)
}


if (exists("snakemake", inherits = FALSE)) run_psap_scoring(workflow_config())
