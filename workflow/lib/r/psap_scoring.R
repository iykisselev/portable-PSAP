# Per-sample lookup selection and PSAP scoring.  Scripts supply only bindings.

SCORING_OUTPUT_COLUMNS <- c("Key", "unit", "partition_id", "Dz.Model", "PSAP", "selected_score_bin", "next_PSAP")


empty_scoring_output <- function() {
  data.table(Key = character(), unit = character(), partition_id = character(), `Dz.Model` = character(), PSAP = numeric(), selected_score_bin = integer(), next_PSAP = numeric())
}


read_scoring_variants <- function(path) {
  variants <- read_tsv(path)
  require_columns(variants, c("Key", "CHROM", "CADD_PHRED", "gt_class", "unit", "partition_id", "partition_status"), "PSAP-ready variant table")
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


lookup_psap <- function(
    variants, lookup, score_max, score_column = "CADD_PHRED",
    lookup_contract = PARTITIONED_NULL_BUNDLE_CONTRACT, metadata = NULL, row_index = NULL) {
  score <- pmin(pmax(variants[[score_column]], 0), score_max)
  if (is_supported_generated_bundle_contract(lookup_contract)) {
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
  setorderv(variants, c(".testing_unit_id", ".partition_id", "CADD_PHRED", "Key"), c(1, 1, -1, 1), na.last = TRUE)
  variants[, .SD[CADD_PHRED == max(CADD_PHRED)], by = .(.testing_unit_id, .partition_id)]
}


score_model <- function(variants, model, genotype_class, lookup, lookup_source) {
  selected <- variants[
    partition_status == "eligible" & gt_class == genotype_class & !is.na(CADD_PHRED) & !is.na(.testing_unit_id) & nzchar(.testing_unit_id) & !is.na(.partition_id) & nzchar(.partition_id)
  ]
  if (nrow(selected) == 0) return(NULL)
  selected <- select_highest(selected)
  row_index <- require_lookup_units(selected, lookup, model)
  selected[, selected_score_bin := generated_lookup_bin(CADD_PHRED, lookup_source$metadata$score_max, lookup_source$metadata$bin_width)]
  selected[, PSAP := lookup_psap(
    selected, lookup, lookup_source$score_max, lookup_contract = lookup_source$type,
    metadata = lookup_source$metadata, row_index = row_index
  )]
  selected[, .(
    Key = as.character(Key), unit = as.character(.testing_unit_id), partition_id = as.character(.partition_id),
    `Dz.Model` = as.character(model), PSAP = as.numeric(PSAP), selected_score_bin = as.integer(selected_score_bin),
    next_PSAP = as.numeric(ifelse(selected_score_bin == ncol(lookup$values), 0, lookup$values[cbind(row_index, pmin(selected_score_bin + 1L, ncol(lookup$values)))]))
  )]
}


variants_for_model <- function(variants, definition) {
  chromosome <- canonical_chromosome(variants$CHROM)
  if (!is.na(definition$required_chromosome)) {
    return(variants[chromosome == definition$required_chromosome])
  }
  variants[chromosome != "X"]
}


score_variants <- function(variants, lookup_source, hemizygote_model) {
  variants <- copy(variants)
  variants[, `:=`(.testing_unit_id = clean_identifier(unit), .partition_id = clean_identifier(partition_id))]
  variants[, gt_class := tolower(trimws(as.character(gt_class)))]
  definitions <- DISEASE_MODELS[hemizygote_model | !DISEASE_MODELS$is_hem, , drop = FALSE]
  results <- list()
  for (current_partition_id in lookup_source$partition_ids) {
    partition_variants <- variants[.partition_id == current_partition_id]
    for (index in seq_len(nrow(definitions))) {
      definition <- definitions[index, ]
      selected <- variants_for_model(partition_variants, definition)
      if (nrow(selected) == 0L) next
      lookup_name <- paste(current_partition_id, definition$model, sep = "\t")
      # Keep only one partition/model matrix resident at a time. Result rows
      # are scalar values, so the matrix can be reclaimed before the next one.
      lookup <- read_lookup(
        lookup_source$paths[[lookup_name]], lookup_source$metadata, lookup_name
      )
      result <- score_model(
        selected, definition$model, definition$genotype_class, lookup, lookup_source
      )
      if (!is.null(result)) results[[length(results) + 1L]] <- result
      rm(lookup)
      gc(verbose = FALSE)
    }
  }
  if (length(results) == 0L) return(empty_scoring_output())
  output <- rbindlist(results, use.names = TRUE, fill = FALSE)
  # Reading by lookup role must not change model/tie order in the public reports.
  output <- output[order(match(output$Dz.Model, definitions$model), output$partition_id)]
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
