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
