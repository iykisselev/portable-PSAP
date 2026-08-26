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
    value[is.na(value) | value %in% c("", ".", "-", "NA", "NaN")] <- NA_character_
    value
}
