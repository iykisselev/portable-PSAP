#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_dir <- normalizePath(file.path(test_dir, ".."))

source(file.path(repo_dir, "workflow", "lib", "r", "common.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "table_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genomic_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "genotypes.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "disease_models.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "bundle_contracts.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "lookup_contracts.R"))
stopifnot(is_supported_generated_bundle_contract("generated_bundle_v3"))
source(file.path(repo_dir, "tests", "r", "null_distribution_oracle.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "psap_scoring.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "candidate_reporting.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "control_reference.R"))
source(file.path(repo_dir, "workflow", "lib", "r", "cohort_reporting.R"))

roundtrip_path <- tempfile(fileext = ".tsv")
write_tsv(data.table(identifier = "one", value = 1.5), roundtrip_path)
roundtrip <- read_tsv(roundtrip_path)
stopifnot(identical(names(roundtrip), c("identifier", "value")), identical(roundtrip$identifier, "one"))
unlink(roundtrip_path)
stopifnot(
    identical(compatible_control_models("DOM-het"), disease_model_names()),
    identical(compatible_control_models("REC-hom"), "REC-hom"),
    identical(standard_disease_models()$model, disease_model_names()[1:4]),
    identical(hem_disease_model()$required_chromosome, "X")
)


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


stopifnot(isTRUE(require_boolean(TRUE, "model.hemizygote")))
expect_error(
    require_boolean("false", "model.hemizygote"),
    "must be a boolean"
)


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

decimal_bins <- fixed_width_bin_index(
    c(0, 0.05, 0.149999, 0.15, 0.150001, 0.3),
    score_max = 0.3,
    bin_width = 0.05
)
stopifnot(
    identical(decimal_bins, c(1L, 2L, 3L, 4L, 4L, 6L)),
    identical(
        generated_bin_index(c(0.05, 0.15, 0.3), 0.3, 0.05),
        c(2L, 4L, 6L)
    ),
    identical(fixed_width_bin_count(0.3, 0.1), 3L)
)

default_width <- 0.05
default_maximum <- 70
default_count <- fixed_width_bin_count(default_maximum, default_width)
internal_boundaries <- vapply(
    seq_len(default_count - 1L),
    function(index) as.numeric(formatC(index * default_width, format = "f", digits = 2L)),
    numeric(1)
)
stopifnot(
    identical(
        fixed_width_bin_index(internal_boundaries - 1e-9, default_maximum, default_width),
        seq_len(default_count - 1L)
    ),
    identical(
        fixed_width_bin_index(internal_boundaries, default_maximum, default_width),
        seq.int(2L, default_count)
    ),
    identical(
        fixed_width_bin_index(internal_boundaries + 1e-9, default_maximum, default_width),
        seq.int(2L, default_count)
    ),
    identical(
        fixed_width_bin_index(default_maximum, default_maximum, default_width),
        default_count
    )
)

af_events <- observed_af_events(c(0, 0.25))
stopifnot(
    identical(af_events$het, c(0, 0.375)),
    identical(af_events$hom, c(0, 0.0625)),
    identical(zero_unobserved_events(2)$het, c(0, 0))
)

canonical_states <- data.table(
    CHROM = c("1", "1", "1"),
    POS = c(1L, 2L, 3L),
    REF = c("A", "A", "A"),
    ALT = c("C", "G", "T"),
    allele_state = c("observed", "absent", "rejected"),
    af = c(0.25, NA_real_, NA_real_),
    state_reason = c(
        "valid_configured_af",
        "not_present_in_population_panel",
        "non_pass_filter"
    )
)
canonical_membership <- data.table(
    testing_unit_id = c("UNIT_A", "UNIT_A", "UNIT_A", "UNIT_B"),
    CHROM = "1",
    POS = c(1L, 2L, 3L, 1L),
    REF = "A",
    ALT = c("C", "G", "T", "C")
)
canonical_scores <- data.table(
    CHROM = "1",
    POS = c(1L, 2L, 3L),
    REF = "A",
    ALT = c("C", "G", "T"),
    score = c(10, 20, 30)
)
canonical_mask <- data.table(CHROM = "1", START = 0L, END = 3L)
assembled_zero <- assemble_calibration_events(
    canonical_states,
    canonical_membership,
    canonical_scores,
    canonical_mask,
    "zero"
)
stopifnot(
    nrow(assembled_zero) == 3L,
    !any(assembled_zero$allele_state == "rejected"),
    all(assembled_zero[allele_state == "absent", effective_af] == 0),
    identical(
        sort(assembled_zero[POS == 1L, testing_unit_id]),
        c("UNIT_A", "UNIT_B")
    )
)
assembled_mutation <- assemble_calibration_events(
    canonical_states,
    canonical_membership,
    canonical_scores,
    canonical_mask,
    "mutation_rate",
    data.table(
        CHROM = c("1", "1"),
        POS = c(2L, 3L),
        REF = "A",
        ALT = c("G", "T"),
        mu = c(1e-8, 0.5)
    )
)
stopifnot(
    nrow(assembled_mutation) == 3L,
    !any(assembled_mutation$POS == 3L),
    identical(
        assembled_mutation[allele_state == "absent", assignment_source],
        "absent_mutation_rate"
    ),
    assembled_mutation[allele_state == "absent", effective_af] == 1e-8
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

y_hem <- apply_sex_classification(
    standardize_genotypes(data.table(CHROM = "Y", GT = "1"), hemizygote_model = TRUE),
    sex = 1,
    hemizygote_model = TRUE
)
stopifnot(nrow(y_hem$genotypes) == 0L)
unknown_x <- apply_sex_classification(
    standardize_genotypes(data.table(CHROM = "X", GT = "0/1"), hemizygote_model = TRUE),
    sex = NA_real_,
    hemizygote_model = TRUE
)
stopifnot(nrow(unknown_x$genotypes) == 0L, unknown_x$removed_unknown_sex_xy == 1L)
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
    list(
        contract = "generated_bundle_v2",
        binning = list(
            algorithm = "fixed_width_decimal_v1",
            score_max = 6,
            bin_width = 2,
            bin_count = 3
        )
    ),
    metadata_path,
    auto_unbox = TRUE
)
validation_path <- tempfile(fileext = ".json")
jsonlite::write_json(
    list(status = "ok", metadata_sha256 = "fixture-digest"),
    validation_path,
    auto_unbox = TRUE
)
marker_path <- tempfile(fileext = ".ok")
writeLines("generated_bundle_v2\tfixture-digest\tfixture-validation-digest", marker_path)
generated_contract <- resolve_lookup_contract(
    "generated_bundle_v2",
    list(het = generated_path),
    metadata_path,
    validation_path,
    marker_path,
    expected_score_max = 6,
    expected_bin_width = 2
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
prepared_lookup_paths <- list(
    het = tempfile(fileext = ".rds"),
    hom = tempfile(fileext = ".rds"),
    hem_chrX = tempfile(fileext = ".rds")
)
prepare_psap_lookups(
    list(
        paths = list(het = generated_path, hom = generated_path, hem_chrX = generated_path),
        metadata = generated_contract$metadata
    ),
    prepared_lookup_paths,
    hemizygote_model = TRUE
)
stopifnot(
    identical(read_prepared_lookup(prepared_lookup_paths$het)$row_of, c(`unit-1` = 1L))
)
hem_scored <- score_variants(
    data.table(
        Key = "X:2:A>C",
        CHROM = "X",
        CADD_PHRED = 2,
        gt_class = "hem",
        unit = "unit-1"
    ),
    list(
        type = "generated_bundle_v2",
        paths = prepared_lookup_paths,
        score_max = 6,
        metadata = list(score_max = 6, bin_width = 2, bin_count = 3, algorithm = "fixed_width_decimal_v1")
    ),
    hemizygote_model = TRUE
)
stopifnot(
    nrow(hem_scored) == 1L,
    identical(hem_scored$Dz.Model, "X-linked-hem"),
    identical(hem_scored$PSAP, 0.75)
)
expect_error(
    resolve_lookup_contract("unknown", list(), NULL),
    "Unsupported lookup contract"
)

candidate_fixture <- build_candidates(
    data.table(
        unit = c("unit-b", "unit-a"),
        model = c("DOM-het", "DOM-het"),
        variant_key = c("1:2:A>C", "1:1:A>G"),
        PSAP = c(0.5, 0.2)
    ),
    "affected", 2
)
stopifnot(
    identical(candidate_fixture$candidate_id, c("affected:DOM:001", "affected:DOM:002")),
    identical(candidate_fixture$variant_key, c("1:1:A>G", "1:2:A>C"))
)

joined_control <- attach_psap(
    data.table(
      pid = "control", phenotype = "unaffected", sex = 1, unit = "unit-a",
      model = "DOM-het", variant_key = "1:1:A>G", PSAP = NA_real_
    ),
    data.table(unit = "unit-a", model = "DOM-het", variant_key = "1:1:A>G", PSAP = 0.25)
)
stopifnot(identical(joined_control$PSAP, 0.25))

report_fixture <- validate_candidates(
    candidate_fixture,
    data.table(
      pid = "control", phenotype = "unaffected", sex = 1, unit = "unit-a",
      model = "DOM-het", variant_key = "1:1:A>G", PSAP = 0.25
    )
)
stopifnot(
    report_fixture$detail[variant_key == "1:1:A>G", validation] == "control_match",
    report_fixture$detail[variant_key == "1:2:A>C", validation] == "ok"
)

compatibility_candidates <- data.table(
    candidate_id = c("dom", "rec-match", "rec-mismatch"),
    pid = "affected",
    phenotype = "affected",
    sex = 2,
    unit = "unit-a",
    model = c("DOM-het", "REC-hom", "REC-hom"),
    variant_key = c("1:3:A>G", "1:3:A>G", "1:4:A>C"),
    PSAP = 0.1
)
compatibility_controls <- data.table(
    pid = c("rec-control", "dom-control"),
    phenotype = "unaffected",
    sex = 1,
    unit = "unit-a",
    model = c("REC-hom", "DOM-het"),
    variant_key = c("1:3:A>G", "1:4:A>C"),
    PSAP = 0.2
)
compatibility_report <- validate_candidates(compatibility_candidates, compatibility_controls)
stopifnot(
    all(compatibility_report$detail[candidate_id %in% c("dom", "rec-match"), validation] == "control_match"),
    compatibility_report$detail[candidate_id == "rec-mismatch", validation] == "ok"
)
no_controls_report <- validate_candidates(compatibility_candidates, compatibility_controls[0])
stopifnot(
    all(no_controls_report$detail$validation == "no_controls"),
    all(no_controls_report$detail$control_count == 0L),
    all(is.na(no_controls_report$detail$control_pids))
)
