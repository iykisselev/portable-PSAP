abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}


require_boolean <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    display_value <- if (length(value)) paste(as.character(value), collapse = ", ") else "NULL"
    abort("%s must be a boolean; received %s", name, display_value)
  }
  value
}


parse_sex <- function(value) {
  normalized <- trimws(as.character(value))
  missing_values <- c("", "NA", "NaN", ".", "-", "-9", "0")

  if (is.na(normalized) || normalized %in% missing_values) {
    return(NA_real_)
  }

  sex <- suppressWarnings(as.numeric(normalized))
  if (is.na(sex)) {
    abort("sex must be numeric or missing; received %s", value)
  }
  sex
}
