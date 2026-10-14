# Per-sample lookup selection and PSAP scoring.  Scripts supply only bindings.

SCORING_OUTPUT_COLUMNS <- c("Key", "unit", "Dz.Model", "PSAP")


empty_scoring_output <- function() {
  data.table(Key = character(), unit = character(), `Dz.Model` = character(), PSAP = numeric())
}


read_scoring_variants <- function(path) {
  variants <- read_tsv(path)
  require_columns(variants, c("Key", "CHROM", "CADD_PHRED", "gt_class", "unit"), "PSAP-ready variant table")
  variants[, CADD_PHRED := suppressWarnings(as.numeric(CADD_PHRED))]
  variants
}


read_lookup <- function(path, metadata = NULL, role = NULL) {
  lookup <- read_tsv(path, header = FALSE, colClasses = "character")
  if (ncol(lookup) < 2 || nrow(lookup) == 0) {
    abort("PSAP lookup table must contain an ID column and score columns: %s", path)
  }
  ids <- trimws(as.character(lookup[[1]]))
  values <- lookup[, -1, with = FALSE]
  value_matrix <- matrix(
    suppressWarnings(as.numeric(unlist(values, use.names = FALSE))),
    nrow = nrow(values), ncol = ncol(values)
  )
  invalid_rows <- rowSums(is.na(value_matrix)) > 0L
  if (identical(tolower(ids[[1L]]), "testing_unit_id") || (invalid_rows[1] && any(!invalid_rows[-1]))) {
    ids <- ids[-1]
    value_matrix <- value_matrix[-1, , drop = FALSE]
  }
  if (anyNA(value_matrix)) abort("PSAP lookup table contains non-numeric score values: %s", path)
  if (anyDuplicated(ids)) abort("PSAP lookup table contains duplicate testing-unit IDs: %s", path)
  lookup <- list(
    ids = ids,
    values = value_matrix,
    row_of = setNames(seq_along(ids), ids)
  )
  if (!is.null(metadata)) {
    validate_generated_lookup(lookup, metadata, if (is.null(role)) "selected" else role)
  }
  lookup
}


prepare_psap_lookups <- function(lookup_source, output_paths, hemizygote_model) {
  roles <- model_roles(hemizygote_model)
  for (role in roles) {
    output <- output_paths[[role]]
    if (is.null(output) || !nzchar(output)) {
      abort("Missing prepared output path for %s lookup", role)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      read_lookup(lookup_source$paths[[role]], lookup_source$metadata, role),
      output,
      version = 2
    )
  }
  invisible(output_paths)
}


read_prepared_lookup <- function(path) {
  lookup <- readRDS(path)
  if (!is.list(lookup) || !all(c("ids", "values", "row_of") %in% names(lookup))) {
    abort("Prepared PSAP lookup has an invalid structure: %s", path)
  }
  lookup
}


lookup_psap <- function(
    variants, lookup, score_max, score_column = "CADD_PHRED",
    lookup_contract = "external_legacy_v0", metadata = NULL, row_index = NULL) {
  score <- pmin(pmax(variants[[score_column]], 0), score_max)
  if (identical(lookup_contract, "external_legacy_v0")) {
    score_grid <- seq(0, score_max, length.out = ncol(lookup$values))
    score_bin <- findInterval(score, score_grid)
    score_bin[score_bin < 1] <- 1
    score_bin[score_bin > ncol(lookup$values)] <- ncol(lookup$values)
  } else if (is_supported_generated_bundle_contract(lookup_contract)) {
    score_bin <- generated_lookup_bin(score, metadata$score_max, metadata$bin_width)
  } else {
    abort("Unsupported lookup contract: %s", lookup_contract)
  }
  if (is.null(row_index)) row_index <- lookup_row_index(lookup, variants$.testing_unit_id)
  if (anyNA(row_index)) {
    missing_ids <- head(unique(variants$.testing_unit_id[is.na(row_index)]), 5)
    abort("PSAP lookup table is missing testing-unit ID(s): %s", paste(missing_ids, collapse = ", "))
  }
  lookup$values[cbind(row_index, score_bin)]
}


select_highest <- function(variants) {
  setorderv(variants, c(".testing_unit_id", "CADD_PHRED", "Key"), c(1, -1, 1), na.last = TRUE)
  variants[, .SD[CADD_PHRED == max(CADD_PHRED)], by = .testing_unit_id]
}


score_model <- function(variants, model, genotype_class, lookup, lookup_source) {
  selected <- variants[
    gt_class == genotype_class & !is.na(CADD_PHRED) & !is.na(.testing_unit_id) & nzchar(.testing_unit_id)
  ]
  if (nrow(selected) == 0) return(NULL)
  selected <- select_highest(selected)
  row_index <- require_lookup_units(selected, lookup, model)
  selected[, PSAP := lookup_psap(
    selected, lookup, lookup_source$score_max, lookup_contract = lookup_source$type,
    metadata = lookup_source$metadata, row_index = row_index
  )]
  selected[, .(
    Key = as.character(Key), unit = as.character(.testing_unit_id),
    `Dz.Model` = as.character(model), PSAP = as.numeric(PSAP)
  )]
}


score_lookup_models <- function(variants, definitions, lookup_source, lookup_name) {
  # This scope owns the only dense lookup. Results contain scalar PSAPs, never
  # a reference to the matrix, so the next role can reclaim this allocation.
  lookup <- read_prepared_lookup(lookup_source$paths[[lookup_name]])
  lapply(seq_len(nrow(definitions)), function(index) {
    definition <- definitions[index, ]
    selected <- variants
    if (definition$is_hem) {
      selected <- variants[canonical_chromosome(CHROM) == definition$required_chromosome]
    }
    score_model(selected, definition$model, definition$genotype_class, lookup, lookup_source)
  })
}


score_variants <- function(variants, lookup_source, hemizygote_model) {
  variants <- copy(variants)
  variants[, .testing_unit_id := clean_identifier(unit)]
  variants[, gt_class := tolower(trimws(as.character(gt_class)))]
  definitions <- DISEASE_MODELS[hemizygote_model | !DISEASE_MODELS$is_hem, , drop = FALSE]
  results <- list()
  for (lookup_name in unique(definitions$lookup)) {
    selected_models <- definitions[definitions$lookup == lookup_name, , drop = FALSE]
    results <- c(results, score_lookup_models(variants, selected_models, lookup_source, lookup_name))
    # Reclaim the previous matrix before allocating another one. Merely replacing
    # a list entry allows unused matrices to accumulate until R chooses to collect.
    gc(verbose = FALSE)
  }
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) return(empty_scoring_output())
  output <- rbindlist(results, use.names = TRUE, fill = FALSE)
  # Reading by lookup role must not change model/tie order in the public reports.
  output <- output[order(match(output$Dz.Model, definitions$model))]
  setcolorder(output, SCORING_OUTPUT_COLUMNS)
  output
}


run_psap_scoring <- function(config) {
  setDTthreads(config$threads)
  variants <- read_scoring_variants(config$input)
  output <- score_variants(variants, lookup_source = config$lookup_source, hemizygote_model = config$hemizygote_model)
  if (nrow(output) > 0) setorderv(output, "PSAP", 1L, na.last = TRUE)
  if (!identical(names(output), SCORING_OUTPUT_COLUMNS)) {
    abort(
      "Internal scoring output schema mismatch: expected %s; found %s",
      paste(SCORING_OUTPUT_COLUMNS, collapse = ", "), paste(names(output), collapse = ", ")
    )
  }
  write_tsv(output, config$output, threads = config$threads)
  message("Candidate rows written: ", nrow(output))
  message("Wrote PSAP results: ", config$output)
}
