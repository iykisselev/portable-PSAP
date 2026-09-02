#!/usr/bin/env Rscript

# Snakemake adapter for unaffected control-reference construction.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) source(file.path(snakemake@scriptdir, "bootstrap.R"))
options(scipen = 999)


workflow_config <- function() {
  inputs <- as.character(snakemake@input[["annotated"]])
  psap_inputs <- as.character(snakemake@input[["psap"]])
  pids <- as.character(snakemake@params[["pids"]])
  sexes <- as.character(snakemake@params[["sexes"]])
  lengths <- c(input = length(inputs), psap_input = length(psap_inputs), pid = length(pids), sex = length(sexes))
  if (length(unique(lengths)) != 1L) abort("Control input counts differ: %s", paste(names(lengths), lengths, sep = "=", collapse = ", "))
  list(
    input = inputs, psap_input = psap_inputs, pid = pids, sex_value = vapply(sexes, parse_sex, numeric(1)),
    hemizygote_model = parse_bool(as.character(snakemake@params[["hemizygote"]]), "model.hemizygote"),
    deduplicate = isTRUE(snakemake@params[["deduplicate"]]), output = snakemake@output[["control_variants"]], threads = snakemake@threads
  )
}


if (exists("snakemake", inherits = FALSE)) run_control_reference(workflow_config())
