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
    memberships = snakemake@input[["memberships"]],
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
  columns <- c("CHROM", "POS", "REF", "ALT", "CADD_PHRED", "GT")
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)

  if (length(readLines(connection, n = 1L, warn = FALSE)) == 0L) {
    return(data.table(
      CHROM = character(), POS = integer(), REF = character(), ALT = character(),
      CADD_PHRED = numeric(), GT = character()
    ))
  }

  genotypes <- read_tsv(path, header = FALSE, col.names = columns)
  if (ncol(genotypes) != length(columns)) {
    abort(
      "Expected six input columns: %s; found %d",
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


expand_units <- function(genotypes, memberships) {
  if (nrow(genotypes) == 0L) {
    genotypes[, `:=`(
      unit = character(), partition_id = character(), partition_status = character(),
      failed_roles_json = character()
    )]
    return(genotypes)
  }
  memberships <- read_tsv(memberships)
  require_columns(memberships, c("CHROM", "POS", "REF", "ALT", "unit", "partition_id", "partition_status", "failed_roles_json"), "Certified partition membership")
  memberships[, `:=`(POS = as.integer(POS), unit = clean_identifier(unit), partition_id = clean_identifier(partition_id))]
  expanded <- merge(genotypes, memberships[, .(CHROM, POS, REF, ALT, unit, partition_id, partition_status, failed_roles_json)],
    by = c("CHROM", "POS", "REF", "ALT"), all = FALSE, sort = FALSE)
  unique(expanded[!is.na(unit) & !is.na(partition_id)], by = c("Key", "unit", "partition_id"))
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

  annotated <- expand_units(classified$genotypes, config$memberships)
  output_columns <- c(
    "Key", "CADD_PHRED", "unit", "partition_id", "partition_status", "failed_roles_json",
    "CHROM", "POS", "REF", "ALT", "GT", "gt_class"
  )
  write_tsv(
    annotated[, output_columns, with = FALSE], config$output,
    threads = config$threads
  )
  message("Annotated rows written: ", nrow(annotated))
  message("Wrote prepared genotypes: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
  run_sample_preparation(workflow_config())
}
