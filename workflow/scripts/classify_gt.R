#!/usr/bin/env Rscript
# Apply sex-aware genotype classification to standardized genotype calls.

suppressPackageStartupMessages({
  library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


# Build the invocation contract and validate the sex/model values before data I/O.
workflow_config <- function() {
  list(
    input = snakemake@input[["gt"]],
    output = snakemake@output[["gt"]],
    sex_value = parse_sex(as.character(snakemake@params[["sex"]])),
    hemizygote = parse_bool(
      as.character(snakemake@params[["hemizygote_model"]]),
      "model.hemizygote"
    ),
    threads = snakemake@threads
  )
}


# Load the standardized genotype table required by the sex-classification helper.
read_genotypes <- function(path) {
  genotypes <- fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    data.table = TRUE,
    showProgress = FALSE
  )

  require_columns(
    genotypes,
    c("Key", "CHROM", "GT", "gt_class"),
    "Standardized genotype table"
  )

  genotypes
}


# Apply the configured sex model and write only the classified genotype records.
run <- function(config) {
  setDTthreads(config$threads)

  message("Reading standardized genotypes: ", config$input)
  genotypes <- read_genotypes(config$input)

  classified <- apply_sex_classification(
    genotypes,
    sex = config$sex_value,
    hemizygote_model = config$hemizygote
  )

  message(
    "Male X/Y heterozygous calls removed: ",
    classified$removed_male_xy_het
  )
  message("Unknown-sex X/Y calls removed: ", classified$removed_unknown_sex_xy)
  message("Output rows: ", nrow(classified$genotypes))

  fwrite(
    classified$genotypes,
    file = config$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    nThread = config$threads
  )

  message("Wrote classified genotypes: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
  run(workflow_config())
}
