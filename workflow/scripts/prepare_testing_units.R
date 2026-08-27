#!/usr/bin/env Rscript

# Normalize the configured cohort-level testing-unit annotations.

suppressPackageStartupMessages({
  library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


# Build the Snakemake input/output contract for this cohort-wide normalization.
workflow_config <- function() {
  list(
    input = snakemake@input[["annotations"]],
    output = snakemake@output[["unit_table"]],
    threads = snakemake@threads
  )
}


# Normalize missing or placeholder testing-unit identifiers to character NA.
clean_unit_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "." | x == "-"] <- NA_character_
  x
}


# Expand comma-delimited unit memberships so every output row has one unit.
expand_units <- function(annotations) {
  if (nrow(annotations) == 0L) {
    annotations[, unit := character()]
    return(annotations)
  }

  pieces <- strsplit(
    clean_unit_value(annotations$unit),
    split = ",",
    fixed = TRUE
  )
  # Repeat source-row positions by their split length to preserve all other
  # annotation fields while creating one row for each unit membership.
  row_index <- rep.int(seq_len(nrow(annotations)), lengths(pieces))
  expanded <- annotations[row_index]
  expanded[, unit := clean_unit_value(unlist(pieces, use.names = FALSE))]
  expanded[!is.na(unit)]
}


# Return one CADD-scored row for each unique variant and testing-unit pairing.
prepare_testing_units <- function(annotations) {
  require_columns(
    annotations,
    c("CHROM", "POS", "REF", "ALT", "CADD_PHRED", "unit"),
    "Configured testing-unit annotation table"
  )

  annotations <- add_vcf_variant_identity(annotations)
  annotations <- expand_units(annotations)
  # The same variant can overlap multiple units, but repeated annotations for
  # one Key/unit pair would duplicate downstream sample evidence.
  annotations <- unique(annotations, by = c("Key", "unit"))
  annotations[, .(Key, CADD_PHRED, unit)]
}


# Normalize the cohort annotations and persist the scoring-table input contract.
run <- function(config) {
  setDTthreads(config$threads)

  message("Reading configured testing-unit annotations: ", config$input)
  annotations <- fread(
    config$input,
    sep = "\t",
    header = TRUE,
    quote = "",
    data.table = TRUE,
    showProgress = FALSE
  )
  input_rows <- nrow(annotations)
  output <- prepare_testing_units(annotations)

  fwrite(
    output,
    file = config$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    nThread = config$threads
  )

  message("Input annotation rows: ", input_rows)
  message("Testing-unit rows written: ", nrow(output))
  message("Unique variants written: ", uniqueN(output$Key))
  message("Unique testing units written: ", uniqueN(output$unit))
  message("Wrote testing-unit table: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
  run(workflow_config())
}
