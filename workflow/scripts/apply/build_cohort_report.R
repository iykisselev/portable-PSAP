#!/usr/bin/env Rscript

# Snakemake adapter for cohort-level candidate/control reporting.

suppressPackageStartupMessages(library(data.table))
if (exists("snakemake", inherits = FALSE)) {
  source(file.path(snakemake@scriptdir, "..", "..", "shared", "bootstrap.R"))
}
options(scipen = 999)


workflow_config <- function() {
  list(
    candidates = as.character(snakemake@input[["candidates"]]),
    evidence = as.character(snakemake@input[["evidence"]]),
    controls = snakemake@input[["controls"]],
    partition_statuses = snakemake@input[["partition_statuses"]],
    output_candidates = snakemake@output[["candidates"]],
    output_summary = snakemake@output[["summary"]],
    output_report_metadata = snakemake@output[["report_metadata"]],
    output_partition_summaries = as.character(snakemake@output[["partition_summaries"]]),
    partition_ids = as.character(snakemake@params[["partition_ids"]]),
    affected_pids = as.character(snakemake@params[["affected_pids"]]),
    affected_sexes = as.character(snakemake@params[["affected_sexes"]]),
    threads = snakemake@threads
  )
}


if (exists("snakemake", inherits = FALSE)) {
  config <- workflow_config()
  candidates <- data.table::rbindlist(lapply(config$candidates, read_tsv), use.names = TRUE, fill = TRUE)
  triaged <- validate_candidates(candidates, read_controls(config$controls))$detail
  reports <- partitioned_cohort_reports(
    triaged, config$partition_ids,
    available_evidence = config$evidence,
    affected_pids = config$affected_pids,
    affected_sexes = config$affected_sexes,
    partition_statuses = config$partition_statuses
  )
  write_tsv(triaged, config$output_candidates, threads = config$threads)
  write_tsv(reports$main, config$output_summary, threads = config$threads)
  for (i in seq_along(config$partition_ids)) write_tsv(reports$partitions[[i]], config$output_partition_summaries[[i]], threads = config$threads)
  jsonlite::write_json(list(
    schema_version = 1L,
    score_direction = "smaller_PSAP_is_stronger",
    partition_ids = as.list(config$partition_ids),
    aggregation_method = "ACAT",
    aggregation_weights = "equal_among_contributing_partitions",
    availability_policy = "compatible_evidence_universe_v1",
    aggregation_input_policy = "raw_survival_v1",
    zero_precedence = "validated_raw_zero_returns_zero_before_endpoint_handling",
    upper_tail_policy = "diagnostic_upper_endpoint_unresolved",
    endpoint_statuses = list(
      zero = c("zero_partition_psap", "zero_and_one_partition_psap"),
      one = c("single_partition_one", "upper_endpoint_unresolved")
    ),
    near_one_thresholds = c(0.99, 0.999999, 0.999999999),
    output_counts = list(main_rows = nrow(reports$main), partition_rows = vapply(reports$partitions, nrow, integer(1)))
  ), config$output_report_metadata, auto_unbox = TRUE, digits = NA, na = "null", pretty = TRUE)
}
