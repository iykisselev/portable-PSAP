#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_dir <- normalizePath(file.path(test_dir, ".."))

source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "acat.R"))

expect_equal <- function(actual, expected, tolerance = 1e-12) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
}

expect_equal(acat_interior(c(0.5, 0.5)), 0.5)
expect_equal(acat_interior(c(0.25, 0.75)), 0.5)
expect_equal(acat_interior(rep(0.2, 3)), 0.2)
expect_equal(acat_interior(c(1e-6, 0.995)), 2.00040004708387e-6, tolerance = 1e-12)

combined <- combine_partition_psaps(c(0, 1, 1e-6))
expect_equal(combined$PSAP, 0)
stopifnot(identical(combined$aggregation_status, "zero_and_one_partition_psap"))
combined <- combine_partition_psaps(1)
expect_equal(combined$PSAP, 1)
stopifnot(identical(combined$aggregation_status, "single_partition_one"))
combined <- combine_partition_psaps(numeric())
stopifnot(is.na(combined$PSAP))
stopifnot(identical(combined$aggregation_status, "no_contributing_partitions"))
combined <- combine_partition_psaps(c(1e-6, 1 - 1e-9))
expect_equal(combined$PSAP, 0.999999997997998, tolerance = 1e-12)
combined <- combine_partition_psaps(c(1e-6, 1))
stopifnot(is.na(combined$PSAP))
stopifnot(identical(combined$aggregation_status, "upper_endpoint_unresolved"))
combined <- combine_partition_psaps(c(1, 1), c(0.995, 0.995), "mid_p")
stopifnot(identical(combined$aggregation_status, "mid_p_exact_one_adjusted"))

message("ACAT tests passed")
