abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}


read_tsv <- function(path, header = TRUE, ...) {
  fread(
    path, sep = "\t", header = header, quote = "", data.table = TRUE,
    showProgress = FALSE, ...
  )
}


write_tsv <- function(data, path, threads = NULL, ...) {
  arguments <- c(
    list(
      x = data, file = path, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = TRUE
    ),
    list(...)
  )
  if (!is.null(threads)) arguments$nThread <- threads
  do.call(fwrite, arguments)
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

  if (is.na(normalized) || tolower(normalized) %in% PED_MISSING_SEX_TOKENS) {
    return(NA_real_)
  }

  sex <- suppressWarnings(as.numeric(normalized))
  if (is.na(sex)) {
    abort("sex must be numeric or missing; received %s", value)
  }
  sex
}
