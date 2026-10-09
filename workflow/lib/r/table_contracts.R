MISSING_TOKENS <- c("", ".", "-", "na", "nan", "none", "null")
PED_MISSING_SEX_TOKENS <- c(MISSING_TOKENS, "-9", "0")


require_columns <- function(x, columns, description) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    abort(
      "%s is missing columns: %s",
      description,
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}


clean_identifier <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value) | tolower(value) %in% MISSING_TOKENS] <- NA_character_
  value
}
