#!/usr/bin/env Rscript

# Construct patient-specific candidates from one affected sample's PSAP output.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))
rm(script_path)

options(scipen = 999)


parse_cli <- function(args) {
    parser <- ArgumentParser(
        description = "Build candidates for one affected sample."
    )

    parser$add_argument("--input", required = TRUE)
    parser$add_argument("--pid", required = TRUE)
    parser$add_argument("--sex", required = TRUE)
    parser$add_argument("--output", required = TRUE)
    parser$add_argument("--threads", type = "integer", default = 1L)

    parsed <- parser$parse_args(args = args)

    require_positive_threads(parsed$threads)

    if (!nzchar(trimws(parsed$pid))) {
        abort("--pid must be non-empty")
    }

    parsed$sex_value <- parse_sex(parsed$sex)
    parsed
}


read_psap <- function(path) {
    psap <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )

    require_columns(
        psap,
        c("Key", "unit", "Dz.Model", "PSAP"),
        "Affected PSAP table"
    )

    psap[, unit := clean_identifier(unit)]
    psap[, model := clean_identifier(`Dz.Model`)]
    psap[, variant_key := clean_identifier(Key)]
    psap[, PSAP := suppressWarnings(as.numeric(PSAP))]

    supported_models <- disease_model_names()
    unsupported <- setdiff(unique(psap$model[!is.na(psap$model)]), supported_models)

    if (length(unsupported) > 0L) {
        abort(
            "Unsupported Dz.Model value(s): %s",
            paste(unsupported, collapse = ", ")
        )
    }

    invalid <- is.na(psap$unit) | is.na(psap$model) | is.na(psap$variant_key)
    if (any(invalid)) {
        abort(
            "Affected PSAP table contains %d row(s) with missing unit, model, or variant key",
            sum(invalid)
        )
    }

    unique(
        psap[, .(unit, model, variant_key, PSAP)],
        by = c("unit", "model", "variant_key", "PSAP")
    )
}


empty_candidates <- function() {
    data.table(
        candidate_id = character(),
        pid = character(),
        phenotype = character(),
        sex = numeric(),
        unit = character(),
        model = character(),
        variant_key = character(),
        PSAP = numeric()
    )
}


build_candidates <- function(psap, pid, sex) {
    eligible <- psap[model %in% disease_model_names()]

    if (nrow(eligible) == 0L) {
        return(empty_candidates())
    }

    setorderv(eligible, c("model", "unit", "variant_key"))
    eligible[, candidate_number := seq_len(.N), by = model]
    eligible[
        ,
        candidate_id := paste(
            pid,
            disease_model_tag(model),
            sprintf("%03d", candidate_number),
            sep = ":"
        )
    ]

    candidates <- eligible[
        ,
        .(
            candidate_id,
            pid = as.character(pid),
            phenotype = "affected",
            sex = sex,
            unit,
            model,
            variant_key,
            PSAP
        )
    ]
    setorderv(
        candidates,
        c("PSAP", "candidate_id", "variant_key"),
        c(1L, 1L, 1L),
        na.last = TRUE
    )
    candidates
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    psap <- read_psap(config$input)
    candidates <- build_candidates(psap, config$pid, config$sex_value)

    fwrite(
        candidates,
        file = config$output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Candidate rows written: ", nrow(candidates))
    message("Wrote affected candidates: ", config$output)
}


if (sys.nframe() == 0) {
    main()
}
