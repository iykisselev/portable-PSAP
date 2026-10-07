#!/usr/bin/env Rscript

# Snakemake adapter for per-sample PSAP scoring.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}
options(scipen = 999)


workflow_config <- function() {
  hemizygote_model <- require_boolean(snakemake@params[["hemizygote"]], "model.hemizygote")
  hem_chrX <- as.character(snakemake@input[["lookup_hem_chrx"]])
  if (length(hem_chrX) == 0L) hem_chrX <- ""
  score_max <- as.numeric(snakemake@params[["score_max"]])
  if (is.na(score_max) || score_max <= 0) abort("model.score_max must be positive")
  if (hemizygote_model && !nzchar(trimws(hem_chrX[[1L]]))) abort("The hemizygous model requires a hem_chrX lookup table")
  lookup_source <- resolve_lookup_contract(
    as.character(snakemake@params[["lookup_contract"]]),
    list(het = snakemake@input[["lookup_het"]], hom = snakemake@input[["lookup_hom"]], hem_chrX = hem_chrX[[1L]]),
    expected_score_max = score_max,
    expected_bin_width = as.character(snakemake@params[["lookup_bin_width"]])
  )
  list(input = snakemake@input[["annotated_gt"]], lookup_source = lookup_source,
       hemizygote_model = hemizygote_model, output = snakemake@output[["psap"]], threads = snakemake@threads)
}


if (exists("snakemake", inherits = FALSE)) run_psap_scoring(workflow_config())
