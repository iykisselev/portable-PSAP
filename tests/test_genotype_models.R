#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_dir <- normalizePath(file.path(test_dir, ".."))

source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genotypes.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "lookup_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "null_distribution.R"))


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

events <- data.table(score = c(10, 20), q = c(0.18, 0.32))
distribution <- conditional_maximum_distribution(events, score_max = 70, bin_width = 5)
stopifnot(
    identical(distribution$status, "ok"),
    isTRUE(all.equal(distribution$event_mass, 0.4424, tolerance = 1e-14)),
    isTRUE(all.equal(distribution$max_mass[c(3, 5)], c(0.1224, 0.32), tolerance = 1e-14)),
    isTRUE(all.equal(distribution$survival[5], 0.32 / 0.4424, tolerance = 1e-14)),
    all(diff(distribution$survival) <= 1e-14)
)

tied <- conditional_maximum_distribution(
    data.table(score = c(10, 20, 20), q = c(0.3, 0.1, 0.2)),
    score_max = 70,
    bin_width = 5
)
stopifnot(
    isTRUE(all.equal(tied$event_mass, 0.496, tolerance = 1e-14)),
    isTRUE(all.equal(tied$max_mass[c(3, 5)], c(0.216, 0.28), tolerance = 1e-14))
)

tiny <- conditional_maximum_distribution(
    data.table(score = 20, q = 1e-24), score_max = 70, bin_width = 5
)
stopifnot(
    identical(tiny$status, "ok"),
    isTRUE(all.equal(tiny$event_mass, 1e-24, tolerance = 1e-35)),
    identical(tiny$survival[5], 1)
)

zero <- conditional_maximum_distribution(
    data.table(score = 20, q = 0), score_max = 70, bin_width = 5
)
stopifnot(identical(zero$status, "zero_event_mass"), zero$event_mass == 0)

af_events <- observed_af_events(c(0, 0.25))
stopifnot(
    identical(af_events$het, c(0, 0.375)),
    identical(af_events$hom, c(0, 0.0625)),
    identical(zero_unobserved_events(2)$het, c(0, 0))
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

source(file.path(repo_dir, "workflow", "scripts", "apply_psap.R"))
tied_scores <- data.table(
    .testing_unit_id = c("unit-1", "unit-1", "unit-1", "unit-2"),
    CADD_PHRED = c(25, 25, 20, 30),
    Key = c("variant-b", "variant-a", "variant-c", "variant-d")
)
selected_scores <- select_highest(copy(tied_scores))
stopifnot(
    identical(
        selected_scores[.testing_unit_id == "unit-1", Key],
        c("variant-a", "variant-b")
    ),
    identical(
        selected_scores[.testing_unit_id == "unit-2", Key],
        "variant-d"
    )
)

generated_path <- tempfile(fileext = ".tsv")
fwrite(
    data.table(testing_unit_id = "unit-1", `1` = 1, `2` = 0.75, `3` = 0.5),
    generated_path,
    sep = "\t"
)
metadata_path <- tempfile(fileext = ".json")
jsonlite::write_json(
    list(contract = "generated_bundle_v1", score_max = 6, bin_width = 2, bin_count = 3),
    metadata_path,
    auto_unbox = TRUE
)
generated_contract <- resolve_lookup_contract(
    "generated_bundle_v1", list(het = generated_path), metadata_path
)
generated_lookup <- read_lookup(generated_path)
generated_scores <- lookup_psap(
    data.table(.testing_unit_id = "unit-1", CADD_PHRED = c(-1, 0, 1.999, 2, 6, 7)),
    generated_lookup,
    score_max = 6,
    lookup_contract = generated_contract$type,
    metadata = generated_contract$metadata
)
stopifnot(identical(generated_scores, c(1, 1, 1, 0.75, 0.5, 0.5)))
expect_error(
    resolve_lookup_contract("unknown", list(), NULL),
    "Unsupported lookup contract"
)
