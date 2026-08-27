#!/usr/bin/env Rscript
# Standardize one-sample genotype calls for downstream PSAP sample preparation.

suppressPackageStartupMessages({
  library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


# Build the invocation contract and validate the configured hemizygote model.
workflow_config <- function() {
  list(
    input = snakemake@input[["gt"]],
    output = snakemake@output[["gt"]],
    hemizygote = parse_bool(
      as.character(snakemake@params[["hemizygote_model"]]),
      "model.hemizygote"
    ),
    threads = snakemake@threads
  )
}


# Read extracted calls, returning a schema-correct empty table for empty input.
read_genotypes <- function(path) {
  columns <- c("CHROM", "POS", "REF", "ALT", "GT")

  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)

  if (length(readLines(connection, n = 1L, warn = FALSE)) == 0L) {
    return(data.table(
      CHROM = character(),
      POS = integer(),
      REF = character(),
      ALT = character(),
      GT = character()
    ))
  }

  genotypes <- fread(
    path,
    sep = "\t",
    header = FALSE,
    col.names = columns,
    quote = "",
    data.table = TRUE,
    showProgress = FALSE
  )

  if (ncol(genotypes) != length(columns)) {
    abort(
      "Expected five input columns: %s; found %d",
      paste(columns, collapse = ", "),
      ncol(genotypes)
    )
  }

  genotypes
}


# Reject conflicting calls and retain one record for each normalized variant key.
deduplicate_by_key <- function(genotypes) {
  # This producer owns the conflicting-call failure and one-row-per-Key contract.
  if (nrow(genotypes) == 0L) {
    return(genotypes)
  }

  calls_by_key <- genotypes[
    ,
    .(
      n_calls = uniqueN(GT),
      calls = paste(sort(unique(GT)), collapse = ",")
    ),
    by = Key
  ]
  conflicts <- calls_by_key[n_calls > 1L]

  if (nrow(conflicts) > 0L) {
    examples <- head(
      paste0(conflicts$Key, " [", conflicts$calls, "]"),
      5L
    )
    abort(
      "Conflicting non-reference calls for %d variant key(s): %s",
      nrow(conflicts),
      paste(examples, collapse = "; ")
    )
  }

  unique(genotypes, by = "Key")
}


# Normalize extracted calls and write the fixed downstream genotype schema.
run <- function(config) {
  setDTthreads(config$threads)

  message("Reading genotype calls: ", config$input)
  genotypes <- read_genotypes(config$input)
  message("Input rows: ", nrow(genotypes))

  genotypes <- add_vcf_variant_identity(genotypes)
  before <- nrow(genotypes)
  genotypes <- standardize_genotypes(
    genotypes,
    hemizygote_model = config$hemizygote
  )
  message("Missing/reference genotype rows removed: ", before - nrow(genotypes))

  before <- nrow(genotypes)
  genotypes <- deduplicate_by_key(genotypes)
  message("Duplicate Key rows removed: ", before - nrow(genotypes))

  output_columns <- c(
    "Key",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "GT",
    "gt_class"
  )

  fwrite(
    genotypes[, output_columns, with = FALSE],
    file = config$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    nThread = config$threads
  )

  message("Output rows: ", nrow(genotypes))
  message("Wrote standardized genotypes: ", config$output)
}


if (exists("snakemake", inherits = FALSE)) {
  run(workflow_config())
}
