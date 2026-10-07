#!/usr/bin/env Rscript

# Parse and validate shared PSAP lookups once before per-sample scoring.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}


workflow_config <- function() {
  hemizygote_model <- require_boolean(snakemake@params[["hemizygote"]], "model.hemizygote")
  score_max <- as.numeric(snakemake@params[["score_max"]])
  if (is.na(score_max) || score_max <= 0) abort("model.score_max must be positive")
  lookup_source <- resolve_lookup_contract(
    as.character(snakemake@params[["lookup_contract"]]),
    list(
      het = snakemake@input[["lookup_het"]],
      hom = snakemake@input[["lookup_hom"]],
      hem_chrX = if (hemizygote_model) snakemake@input[["lookup_hem_chrx"]][[1L]] else ""
    ),
    expected_score_max = score_max,
    expected_bin_width = as.character(snakemake@params[["lookup_bin_width"]])
  )
  list(
    lookup_source = lookup_source,
    output_paths = list(
      het = snakemake@output[["lookup_het"]],
      hom = snakemake@output[["lookup_hom"]],
      hem_chrX = if (hemizygote_model) snakemake@output[["lookup_hem_chrx"]][[1L]] else ""
    ),
    hemizygote_model = hemizygote_model
  )
}


if (exists("snakemake", inherits = FALSE)) {
  config <- workflow_config()
  prepare_psap_lookups(config$lookup_source, config$output_paths, config$hemizygote_model)
}
