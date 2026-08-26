#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_dir <- normalizePath(file.path(test_dir, ".."))

source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genotypes.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "lookup_contracts.R"))


expect_error <- function(expression, pattern) {
    message <- tryCatch(
        {
            force(expression)
            NA_character_
        },
        error = function(condition) conditionMessage(condition)
    )
    stopifnot(!is.na(message), grepl(pattern, message))
}


haploid <- data.table(CHROM = "chrX", GT = "1")
standardized <- standardize_genotypes(copy(haploid), hemizygote_model = TRUE)
stopifnot(
    identical(standardized$GT, "1"),
    identical(standardized$gt_class, "haploid")
)

classified <- apply_sex_classification(
    copy(standardized),
    sex = 1,
    hemizygote_model = TRUE
)
stopifnot(
    nrow(classified$genotypes) == 1L,
    identical(classified$genotypes$gt_class, "hem")
)

reference_haploid <- standardize_genotypes(
    data.table(CHROM = "X", GT = "0"),
    hemizygote_model = TRUE
)
stopifnot(nrow(reference_haploid) == 0L)

expect_error(
    standardize_genotypes(copy(haploid), hemizygote_model = FALSE),
    "hemizygote model is disabled"
)
expect_error(
    standardize_genotypes(
        data.table(CHROM = "1", GT = "1"),
        hemizygote_model = TRUE
    ),
    "only on chromosomes X and Y"
)
expect_error(
    apply_sex_classification(
        copy(standardized),
        sex = 2,
        hemizygote_model = TRUE
    ),
    "female sample"
)
expect_error(
    apply_sex_classification(
        copy(standardized),
        sex = NA_real_,
        hemizygote_model = TRUE
    ),
    "known male sample sex"
)

diploid <- standardize_genotypes(
    data.table(CHROM = c("X", "X"), GT = c("1/1", "0/1")),
    hemizygote_model = TRUE
)
classified_diploid <- apply_sex_classification(
    diploid,
    sex = 1,
    hemizygote_model = TRUE
)
stopifnot(
    nrow(classified_diploid$genotypes) == 1L,
    identical(classified_diploid$genotypes$gt_class, "hem"),
    classified_diploid$removed_male_xy_het == 1L
)

expect_error(
    require_lookup_units(
        data.table(.testing_unit_id = c("unit-1", "missing-unit")),
        list(ids = "unit-1"),
        "DOM-het"
    ),
    "DOM-het lookup is missing 1 observed testing-unit ID"
)
