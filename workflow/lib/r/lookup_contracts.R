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


generated_lookup_bin <- function(score, score_max, bin_width) {
  if (!is.finite(score_max) || score_max <= 0 || !is.finite(bin_width) || bin_width <= 0) {
    abort("Generated lookup metadata requires positive score_max and bin_width")
  }
  count <- ceiling(score_max / bin_width)
  score <- pmin(pmax(as.numeric(score), 0), score_max)
  pmin(floor(score / bin_width) + 1L, count)
}


resolve_lookup_contract <- function(contract, lookup_paths, metadata_path = NULL) {
  if (identical(contract, "external_legacy_v0")) {
    return(list(type = contract, paths = lookup_paths, metadata = NULL))
  }
  if (!identical(contract, "generated_bundle_v1")) {
    abort("Unsupported lookup contract: %s", contract)
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    abort("generated_bundle_v1 requires a metadata.json path")
  }
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = TRUE)
  if (!identical(metadata$contract, "generated_bundle_v1")) {
    abort("Unsupported generated bundle contract: %s", metadata$contract)
  }
  required <- c("score_max", "bin_width", "bin_count")
  if (!all(required %in% names(metadata))) {
    abort("Generated bundle metadata is missing grid fields")
  }
  if (metadata$bin_count != ceiling(metadata$score_max / metadata$bin_width)) {
    abort("Generated bundle metadata has inconsistent bin_count")
  }
  list(type = contract, paths = lookup_paths, metadata = metadata)
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
