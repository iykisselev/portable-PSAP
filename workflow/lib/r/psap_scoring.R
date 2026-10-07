# Per-sample lookup selection and PSAP scoring.  Scripts supply only bindings.

SCORING_OUTPUT_COLUMNS <- c("Key", "unit", "Dz.Model", "PSAP")


empty_scoring_output <- function() {
  data.table(Key = character(), unit = character(), `Dz.Model` = character(), PSAP = numeric())
}


read_scoring_variants <- function(path) {
  variants <- fread(path, sep = "\t", header = TRUE, quote = "", data.table = TRUE, showProgress = FALSE)
  require_columns(variants, c("Key", "CHROM", "CADD_PHRED", "gt_class", "unit"), "PSAP-ready variant table")
  variants[, CADD_PHRED := suppressWarnings(as.numeric(CADD_PHRED))]
  variants
}


read_lookup <- function(path, metadata = NULL, role = NULL) {
  lookup <- fread(path, sep = "\t", header = FALSE, quote = "", colClasses = "character", data.table = TRUE, showProgress = FALSE)
  if (ncol(lookup) < 2 || nrow(lookup) == 0) {
    abort("PSAP lookup table must contain an ID column and score columns: %s", path)
  }
  ids <- trimws(as.character(lookup[[1]]))
  values <- lookup[, -1, with = FALSE]
  value_matrix <- matrix(suppressWarnings(as.numeric(unlist(values, use.names = FALSE))), nrow = nrow(values), ncol = ncol(values))
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
  roles <- c("het", "hom")
  if (hemizygote_model) roles <- c(roles, "hem_chrX")
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


lookup_psap <- function(variants, lookup, score_max, score_column = "CADD_PHRED", lookup_contract = "external_legacy_v0", metadata = NULL, lookup_role = "selected", row_index = NULL) {
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


score_model <- function(variants, model, genotype_class, lookup, lookup_source, lookup_role) {
  selected <- variants[gt_class == genotype_class & !is.na(CADD_PHRED) & !is.na(.testing_unit_id) & nzchar(.testing_unit_id)]
  if (nrow(selected) == 0) return(NULL)
  selected <- select_highest(selected)
  row_index <- require_lookup_units(selected, lookup, model)
  selected[, PSAP := lookup_psap(selected, lookup, lookup_source$score_max, lookup_contract = lookup_source$type, metadata = lookup_source$metadata, lookup_role = lookup_role, row_index = row_index)]
  selected[, .(Key = as.character(Key), unit = as.character(.testing_unit_id), `Dz.Model` = as.character(model), PSAP = as.numeric(PSAP))]
}


score_variants <- function(variants, lookup_source, hemizygote_model) {
  variants <- copy(variants)
  variants[, .testing_unit_id := clean_identifier(unit)]
  variants[, gt_class := tolower(trimws(as.character(gt_class)))]
  definitions <- DISEASE_MODELS[DISEASE_MODELS$lookup != "hem_chrX", ]
  lookup_names <- unique(definitions$lookup)
  lookups <- setNames(
    lapply(lookup_names, function(name) {
      read_prepared_lookup(lookup_source$paths[[name]])
    }),
    lookup_names
  )
  results <- lapply(seq_len(nrow(definitions)), function(index) {
    definition <- definitions[index, ]
    score_model(variants, definition$model, definition$genotype_class, lookups[[definition$lookup]], lookup_source, definition$lookup)
  })
  if (hemizygote_model) {
    hem_variants <- copy(variants)
    hem_variants[, .chromosome := canonical_chromosome(CHROM)]
    hem_definition <- DISEASE_MODELS[DISEASE_MODELS$lookup == "hem_chrX", , drop = FALSE]
    hem_lookup <- read_prepared_lookup(lookup_source$paths[[hem_definition$lookup]])
    results <- c(results, list(score_model(hem_variants[.chromosome == "X" & gt_class == "hem"], hem_definition$model, hem_definition$genotype_class, hem_lookup, lookup_source, hem_definition$lookup)))
  }
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(empty_scoring_output())
  output <- rbindlist(results, use.names = TRUE, fill = FALSE)
  setcolorder(output, SCORING_OUTPUT_COLUMNS)
  output
}


run_psap_scoring <- function(config) {
  setDTthreads(config$threads)
  variants <- read_scoring_variants(config$input)
  output <- score_variants(variants, lookup_source = config$lookup_source, hemizygote_model = config$hemizygote_model)
  if (nrow(output) > 0) setorderv(output, "PSAP", 1L, na.last = TRUE)
  if (!identical(names(output), SCORING_OUTPUT_COLUMNS)) {
    abort("Internal scoring output schema mismatch: expected %s; found %s", paste(SCORING_OUTPUT_COLUMNS, collapse = ", "), paste(names(output), collapse = ", "))
  }
  fwrite(output, file = config$output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, nThread = config$threads)
  message("Candidate rows written: ", nrow(output))
  message("Wrote PSAP results: ", config$output)
}
