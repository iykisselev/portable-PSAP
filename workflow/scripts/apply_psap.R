#!/usr/bin/env Rscript

# Snakemake adapter for per-sample PSAP scoring.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) source(file.path(snakemake@scriptdir, "bootstrap.R"))
options(scipen = 999)


workflow_config <- function() {
  optional_input <- function(name) {
    values <- as.character(snakemake@input[[name]])
    if (length(values)) values else NULL
  }
  hemizygote_model <- parse_bool(as.character(snakemake@params[["hemizygote"]]), "model.hemizygote")
  hem_chrX <- as.character(snakemake@input[["lookup_hem_chrx"]])
  if (length(hem_chrX) == 0L) hem_chrX <- ""
  score_max <- as.numeric(snakemake@params[["score_max"]])
  if (is.na(score_max) || score_max <= 0) abort("model.score_max must be positive")
  if (hemizygote_model && !nzchar(trimws(hem_chrX[[1L]]))) abort("The hemizygous model requires a hem_chrX lookup table")
  lookup_source <- resolve_lookup_contract(
    as.character(snakemake@params[["lookup_contract"]]),
    list(het = snakemake@input[["lookup_het"]], hom = snakemake@input[["lookup_hom"]], hem_chrX = hem_chrX[[1L]]),
    metadata_path = optional_input("lookup_metadata"), validation_path = optional_input("lookup_validation"),
    marker_path = optional_input("lookup_marker"), expected_score_max = score_max,
    expected_bin_width = as.character(snakemake@params[["lookup_bin_width"]]),
    eligible_units = optional_input("eligible_units"), calibration_mask = optional_input("calibration_mask"),
    excluded_alleles = optional_input("excluded_alleles"), excluded_units = optional_input("excluded_units")
  )
  list(input = snakemake@input[["annotated_gt"]], lookup_source = lookup_source,
       hemizygote_model = hemizygote_model, output = snakemake@output[["psap"]], threads = snakemake@threads)
}


if (exists("snakemake", inherits = FALSE)) run_psap_scoring(workflow_config())
