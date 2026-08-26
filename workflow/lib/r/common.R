`%||%` <- function(value, replacement) {
    if (is.null(value)) replacement else value
}


abort <- function(...) {
    stop(sprintf(...), call. = FALSE)
}


parse_bool <- function(value, name) {
    normalized <- tolower(trimws(as.character(value)))
    if (identical(normalized, "true")) return(TRUE)
    if (identical(normalized, "false")) return(FALSE)
    abort("%s must be True or False; received %s", name, value)
}


parse_sex <- function(value) {
    normalized <- trimws(as.character(value))
    missing_values <- c("", "NA", "NaN", ".", "-", "-9", "0")

    if (is.na(normalized) || normalized %in% missing_values) {
        return(NA_real_)
    }

    sex <- suppressWarnings(as.numeric(normalized))
    if (is.na(sex)) {
        abort("--sex must be numeric or missing; received %s", value)
    }
    sex
}


require_positive_threads <- function(value) {
    if (is.na(value) || value < 1L) {
        abort("--threads must be a positive integer")
    }
    invisible(TRUE)
}
