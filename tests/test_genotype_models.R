#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
test_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
repo_dir <- normalizePath(file.path(test_dir, ".."))
source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "model_catalog.R"))
PSAP_MODEL_CATALOG <- read_model_catalog(file.path(repo_dir, "workflow", "model_catalog.json"))
source(file.path(repo_dir, "workflow", "lib", "r", "disease_models.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genomic_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genotypes.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "bundle_contracts.R"))

stopifnot(is_supported_generated_bundle_contract("partitioned_null_bundle_v1"))
stopifnot(identical(disease_model_names(), c("DOM-het", "REC-hom", "X-linked-hom", "X-linked-het", "X-linked-hem")))
male_x <- apply_sex_classification(data.table::data.table(CHROM = "X", gt_class = "hom"), 1, TRUE)$genotypes
stopifnot(identical(male_x$gt_class, "hem"))
message("Partitioned genotype-model tests passed")
