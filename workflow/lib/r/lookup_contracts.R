require_lookup_units <- function(variants, lookup, model) {
  lookup_row <- match(variants$.testing_unit_id, lookup$ids)
  missing <- is.na(lookup_row)

  if (any(missing)) {
    missing_ids <- sort(unique(variants$.testing_unit_id[missing]))
    examples <- paste(head(missing_ids, 10L), collapse = ", ")
    remainder <- length(missing_ids) - min(length(missing_ids), 10L)
    suffix <- if (remainder > 0L) {
      sprintf(" (and %d more)", remainder)
    } else {
      ""
    }
    abort(
      paste0(
        "Model %s lookup is missing %d observed testing-unit ID(s) ",
        "across %d variant row(s): %s%s"
      ),
      model,
      length(missing_ids),
      sum(missing),
      examples,
      suffix
    )
  }

  invisible(TRUE)
}


fixed_width_bin_count <- function(score_max, bin_width) {
  if (!is.finite(score_max) || score_max <= 0 || !is.finite(bin_width) || bin_width <= 0) {
    abort("Generated lookup metadata requires positive score_max and bin_width")
  }
  scaled <- score_max / bin_width
  nearest <- round(scaled)
  tolerance <- 64 * .Machine$double.eps * max(1, abs(scaled))
  if (abs(scaled - nearest) <= tolerance) scaled <- nearest
  as.integer(ceiling(scaled))
}


fixed_width_bin_index <- function(score, score_max, bin_width) {
  count <- fixed_width_bin_count(score_max, bin_width)
  score <- pmin(pmax(as.numeric(score), 0), score_max)
  scaled <- score / bin_width
  nearest <- round(scaled)
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(scaled))
  on_boundary <- abs(scaled - nearest) <= tolerance
  scaled[on_boundary] <- nearest[on_boundary]
  as.integer(pmin(floor(scaled) + 1L, count))
}


generated_lookup_bin <- fixed_width_bin_index


resolve_lookup_contract <- function(
    contract,
    lookup_paths,
    metadata_path = NULL,
    validation_path = NULL,
    marker_path = NULL,
    expected_score_max = NULL,
    expected_bin_width = NULL,
    eligible_units = NULL,
    calibration_mask = NULL,
    excluded_alleles = NULL,
    excluded_units = NULL
) {
  if (identical(contract, "external_legacy_v0")) {
    return(list(
      type = contract,
      paths = lookup_paths,
      metadata = NULL,
      score_max = as.numeric(expected_score_max),
      bin_width = NULL,
      eligible_units = eligible_units,
      calibration_mask = calibration_mask,
      excluded_alleles = excluded_alleles,
      excluded_units = excluded_units
    ))
  }
  if (!is_supported_generated_bundle_contract(contract)) {
    abort("Unsupported lookup contract: %s", contract)
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    abort("Generated bundle requires a metadata.json path")
  }
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = TRUE)
  if (!identical(metadata$contract, contract)) {
    abort("Unsupported generated bundle contract: %s", metadata$contract)
  }
  required <- c("score_max", "bin_width", "bin_count", "algorithm")
  if (is.null(metadata$binning) || !all(required %in% names(metadata$binning))) {
    abort("Generated bundle metadata is missing grid fields")
  }
  grid <- metadata$binning
  if (!identical(grid$algorithm, FIXED_WIDTH_BINNING_ALGORITHM)) {
    abort("Generated bundle uses an unsupported binning algorithm: %s", grid$algorithm)
  }
  if (grid$bin_count != fixed_width_bin_count(grid$score_max, grid$bin_width)) {
    abort("Generated bundle metadata has inconsistent bin_count")
  }
  if (!is.null(expected_score_max) && abs(grid$score_max - as.numeric(expected_score_max)) > 1e-12) {
    abort("Generated bundle score_max differs from the resolved lookup source")
  }
  if (!is.null(expected_bin_width) && nzchar(as.character(expected_bin_width)) &&
      abs(grid$bin_width - as.numeric(expected_bin_width)) > 1e-15) {
    abort("Generated bundle bin_width differs from the resolved lookup source")
  }
  if (is.null(validation_path) || is.null(marker_path) ||
      !file.exists(validation_path) || !file.exists(marker_path)) {
    abort("Generated bundle requires validation.json and validation.ok")
  }
  validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
  marker <- strsplit(trimws(readLines(marker_path, warn = FALSE)), "\t", fixed = TRUE)[[1L]]
  if (!identical(validation$status, "ok") || length(marker) != 3L ||
      marker[[1L]] != contract || marker[[2L]] != validation$metadata_sha256) {
    abort("Generated bundle validation marker is inconsistent")
  }
  normalized_metadata <- list(
    score_max = as.numeric(grid$score_max),
    bin_width = as.numeric(grid$bin_width),
    bin_count = as.integer(grid$bin_count),
    algorithm = grid$algorithm
  )
  list(
    type = contract,
    paths = lookup_paths,
    metadata = normalized_metadata,
    score_max = normalized_metadata$score_max,
    bin_width = normalized_metadata$bin_width,
    eligible_units = eligible_units,
    calibration_mask = calibration_mask,
    excluded_alleles = excluded_alleles,
    excluded_units = excluded_units
  )
}


validate_generated_lookup <- function(lookup, metadata, role) {
  if (ncol(lookup$values) != metadata$bin_count) {
    abort("Generated %s lookup has %d bins; metadata requires %d", role, ncol(lookup$values), metadata$bin_count)
  }
  if (any(!is.finite(lookup$values)) || any(lookup$values < 0 | lookup$values > 1)) {
    abort("Generated %s lookup has values outside [0, 1]", role)
  }
  if (any(apply(lookup$values, 1, function(x) any(diff(x) > 1e-12)))) {
    abort("Generated %s lookup survival values are not nonincreasing", role)
  }
  invisible(TRUE)
}
