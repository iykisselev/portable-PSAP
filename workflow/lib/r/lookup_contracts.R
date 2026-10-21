lookup_row_index <- function(lookup, ids) {
  if (!is.null(lookup$row_of)) {
    return(unname(lookup$row_of[as.character(ids)]))
  }
  match(ids, lookup$ids)
}


require_lookup_units <- function(variants, lookup, model) {
  lookup_row <- lookup_row_index(lookup, variants$.testing_unit_id)
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

  invisible(lookup_row)
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
    expected_score_max = NULL,
    expected_bin_width = NULL
) {
  if (!is_supported_generated_bundle_contract(contract)) {
    abort("Unsupported lookup contract: %s", contract)
  }
  # The validate_model workflow barrier verifies the generated-bundle
  # certificate and its components.  Scoring only needs the already-resolved
  # numeric grid, so do not re-read the certificate for every sample.
  if (is.null(expected_score_max) || !is.finite(as.numeric(expected_score_max)) ||
      as.numeric(expected_score_max) <= 0) {
    abort("Generated lookup source requires a positive resolved score_max")
  }
  if (is.null(expected_bin_width) || !nzchar(as.character(expected_bin_width)) ||
      !is.finite(as.numeric(expected_bin_width)) || as.numeric(expected_bin_width) <= 0) {
    abort("Generated lookup source requires a positive resolved bin_width")
  }
  normalized_metadata <- list(
    score_max = as.numeric(expected_score_max),
    bin_width = as.numeric(expected_bin_width),
    bin_count = fixed_width_bin_count(as.numeric(expected_score_max), as.numeric(expected_bin_width)),
    algorithm = FIXED_WIDTH_BINNING_ALGORITHM
  )
  list(
    type = contract,
    paths = lookup_paths,
    metadata = normalized_metadata,
    score_max = normalized_metadata$score_max,
    bin_width = normalized_metadata$bin_width
  )
}


validate_generated_lookup <- function(lookup, metadata, role) {
  if (ncol(lookup$values) != metadata$bin_count) {
    abort("Generated %s lookup has %d bins; metadata requires %d", role, ncol(lookup$values), metadata$bin_count)
  }
  if (any(!is.finite(lookup$values)) || any(lookup$values < 0 | lookup$values > 1)) {
    abort("Generated %s lookup has values outside [0, 1]", role)
  }
  if (ncol(lookup$values) > 1L &&
      any(lookup$values[, -1L, drop = FALSE] - lookup$values[, -ncol(lookup$values), drop = FALSE] > 1e-12)) {
    abort("Generated %s lookup survival values are not nonincreasing", role)
  }
  invisible(TRUE)
}
