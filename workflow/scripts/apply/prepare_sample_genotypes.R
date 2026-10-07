#!/usr/bin/env Rscript
# Prepare one sample's annotated genotype table without materializing R intermediates.

suppressPackageStartupMessages({
  library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}

options(scipen = 999)


workflow_config <- function() {
  list(
    gt = snakemake@input[["gt"]],
    output = snakemake@output[["annotated_gt"]],
    sex_value = parse_sex(as.character(snakemake@params[["sex"]])),
    hemizygote = require_boolean(
      snakemake@params[["hemizygote_model"]],
      "model.hemizygote"
    ),
    threads = snakemake@threads
  )
}


read_extracted_genotypes <- function(path) {
  columns <- c("CHROM", "POS", "REF", "ALT", "CADD_PHRED", "unit", "GT")
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)

  if (length(readLines(connection, n = 1L, warn = FALSE)) == 0L) {
    return(data.table(
      CHROM = character(), POS = integer(), REF = character(), ALT = character(),
      CADD_PHRED = numeric(), unit = character(), GT = character()
    ))
  }

  genotypes <- fread(
    path, sep = "\t", header = FALSE, col.names = columns, quote = "",
    data.table = TRUE, showProgress = FALSE
  )
  if (ncol(genotypes) != length(columns)) {
    abort(
      "Expected seven input columns: %s; found %d",
      paste(columns, collapse = ", "), ncol(genotypes)
    )
  }
  genotypes
}


deduplicate_by_key <- function(genotypes) {
  if (nrow(genotypes) == 0L) return(genotypes)
  calls_by_key <- genotypes[, .(
    n_calls = uniqueN(GT), calls = paste(sort(unique(GT)), collapse = ",")
  ), by = Key]
  conflicts <- calls_by_key[n_calls > 1L]
  if (nrow(conflicts) > 0L) {
    examples <- head(paste0(conflicts$Key, " [", conflicts$calls, "]"), 5L)
    abort(
      "Conflicting non-reference calls for %d variant key(s): %s",
      nrow(conflicts), paste(examples, collapse = "; ")
    )
  }
  unique(genotypes, by = "Key")
}


clean_unit_value <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value) | value == "" | value == "." | value == "-"] <- NA_character_
  value
}


expand_units <- function(genotypes) {
  if (nrow(genotypes) == 0L) {
    genotypes[, unit := character()]
    return(genotypes)
  }
  pieces <- strsplit(clean_unit_value(genotypes$unit), split = ",", fixed = TRUE)
  expanded <- genotypes[rep.int(seq_len(nrow(genotypes)), lengths(pieces))]
  expanded[, unit := clean_unit_value(unlist(pieces, use.names = FALSE))]
  unique(expanded[!is.na(unit)], by = c("Key", "unit"))
}


run_sample_preparation <- function(config) {
  setDTthreads(config$threads)
  genotypes <- read_extracted_genotypes(config$gt)
  message("Extracted genotype rows: ", nrow(genotypes))

  genotypes <- add_vcf_variant_identity(genotypes)
  before <- nrow(genotypes)
  genotypes <- standardize_genotypes(genotypes, config$hemizygote)
  message("Missing/reference genotype rows removed: ", before - nrow(genotypes))

  before <- nrow(genotypes)
  genotypes <- deduplicate_by_key(genotypes)
  message("Duplicate Key rows removed: ", before - nrow(genotypes))

  classified <- apply_sex_classification(
    genotypes, sex = config$sex_value, hemizygote_model = config$hemizygote
  )
  message("Male X/Y heterozygous calls removed: ", classified$removed_male_xy_het)
  message("Unknown-sex X/Y calls removed: ", classified$removed_unknown_sex_xy)

  annotated <- expand_units(classified$genotypes)
  output_columns <- c(
    "Key", "CADD_PHRED", "unit", "CHROM", "POS", "REF", "ALT", "GT", "gt_class"
  )
  fwrite(
    annotated[, output_columns, with = FALSE], file = config$output, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, nThread = config$threads
  )
  message("Annotated rows written: ", nrow(annotated))
  message("Wrote prepared genotypes: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
  run_sample_preparation(workflow_config())
}
