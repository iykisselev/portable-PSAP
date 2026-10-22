#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
test_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
repo_dir <- normalizePath(file.path(test_dir, ".."))
source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "model_catalog.R"))
PSAP_MODEL_CATALOG <- read_model_catalog(file.path(repo_dir, "workflow", "model_catalog.json"))
source(file.path(repo_dir, "workflow", "lib", "r", "acat.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "table_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "disease_models.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "cohort_reporting.R"))
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

candidates <- data.table::data.table(
  candidate_id = c("P1:DOM:001", "P1:DOM:002"), pid = "P1", phenotype = "affected", sex = 2,
  unit = "G", partition_id = c("cds", "intronic"), model = "DOM-het",
  variant_key = c("1:10_A>T", "1:20_A>G"), PSAP = c(0.02, 0.3),
  selected_score_bin = c(400L, 600L), next_PSAP = c(0.01, 0.2),
  validation = "no_controls", control_count = 0L, control_pids = NA_character_
)
reports <- partitioned_cohort_reports(candidates, c("cds", "intronic"))
stopifnot(nrow(reports$main) == 1L, reports$main$contributing_partition_count == 2L)
stopifnot(identical(reports$main$aggregation_status, "ok"))
stopifnot(all(reports$detail$status == "scored"))
message("Partition reporting tests passed")
