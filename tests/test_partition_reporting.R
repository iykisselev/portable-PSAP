#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
test_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
repo_dir <- normalizePath(file.path(test_dir, ".."))
source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "acat.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "table_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "partition_reporting.R"))

scored <- data.table::data.table(
  pid = c("P1", "P1", "P1"), unit = "G", model = "DOM-het",
  partition_id = c("cds", "cds", "intronic"), PSAP = c(0.1, 0.02, 0.3)
)
partitions <- partition_best_psaps(scored)
stopifnot(partitions[partition_id == "cds", PSAP] == 0.02)
omnibus <- omnibus_partition_psaps(partitions, c("cds", "intronic"))
stopifnot(omnibus$contributing_partition_count == 2L)
stopifnot(identical(omnibus$aggregation_status, "ok"))
stopifnot(grepl('^\\["cds","intronic"\\]$', omnibus$contributing_partitions_json))
message("Partition reporting tests passed")
