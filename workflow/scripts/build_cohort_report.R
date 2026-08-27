#!/usr/bin/env Rscript

# Validate affected candidates against unaffected carried-variant evidence.

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
        description = "Build detailed and summarized PSAP cohort reports."
    )

    parser$add_argument("--candidates", nargs = "+", required = TRUE)
    parser$add_argument("--controls", required = TRUE)
    parser$add_argument("--output-candidates", required = TRUE)
    parser$add_argument("--output-summary", required = TRUE)
    parser$add_argument("--threads", type = "integer", default = 1L)

    parsed <- parser$parse_args(args = args)
    require_positive_threads(parsed$threads)

    parsed
}


read_candidate_tables <- function(paths) {
    required <- c(
        "candidate_id",
        "pid",
        "phenotype",
        "sex",
        "unit",
        "model",
        "variant_key",
        "PSAP"
    )

    tables <- lapply(paths, function(path) {
        candidates <- fread(
            path,
            sep = "\t",
            header = TRUE,
            quote = "",
            data.table = TRUE,
            showProgress = FALSE
        )
        require_columns(candidates, required, paste("Candidate table", path))
        candidates[, required, with = FALSE]
    })

    candidates <- rbindlist(tables, use.names = TRUE, fill = FALSE)
    candidates[, PSAP := suppressWarnings(as.numeric(PSAP))]

    if (any(candidates$phenotype != "affected")) {
        abort("Candidate tables contain rows not labelled affected")
    }

    candidates
}


read_controls <- function(path) {
    required <- c(
        "pid",
        "phenotype",
        "sex",
        "unit",
        "model",
        "variant_key",
        "PSAP"
    )
    controls <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
    )
    require_columns(controls, required, "Control reference")
    controls <- controls[, required, with = FALSE]

    if (any(controls$phenotype != "unaffected")) {
        abort("Control reference contains rows not labelled unaffected")
    }

    unique(controls, by = c("pid", "unit", "model", "variant_key"))
}


minimum_or_na <- function(value) {
    if (length(value) == 0L || all(is.na(value))) return(NA_real_)
    min(value, na.rm = TRUE)
}


collapse_ids <- function(value) {
    value <- unique(value[!is.na(value) & nzchar(value)])
    if (length(value) == 0L) return(NA_character_)
    paste(sort(value), collapse = ",")
}


combine_id_lists <- function(value) {
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value) == 0L) return(NA_character_)
    collapse_ids(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
}


validation_priority <- function(value) {
    if ("control_match" %in% value) return("control_match")
    if (all(value == "no_controls")) return("no_controls")
    "ok"
}


build_candidate_groups <- function(candidates) {
    inconsistent <- candidates[
        ,
        .(
            pid_count = uniqueN(pid),
            unit_count = uniqueN(unit),
            model_count = uniqueN(model)
        ),
        by = candidate_id
    ][pid_count != 1L | unit_count != 1L | model_count != 1L]

    if (nrow(inconsistent) > 0L) {
        abort("candidate_id values span multiple patients, units, or models")
    }

    groups <- candidates[
        ,
        .(
            pid = pid[[1]],
            unit = unit[[1]],
            model = model[[1]],
            variant_keys = list(sort(unique(variant_key))),
            PSAP = minimum_or_na(PSAP)
        ),
        by = candidate_id
    ]
    groups[, variant_count := lengths(variant_keys)]

    if (any(groups$variant_count != 1L)) {
        abort("Candidates must contain exactly one distinct variant")
    }

    groups[
        ,
        candidate_signature := vapply(
            variant_keys,
            paste,
            collapse = "|",
            FUN.VALUE = character(1)
        )
    ]
    groups
}


validate_one_candidate <- function(group, controls) {
    if (nrow(controls) == 0L) {
        return(list(
            validation = "no_controls",
            control_count = 0L,
            control_pids = NA_character_
        ))
    }

    candidate_model <- group$model[[1]]
    candidate_unit <- group$unit[[1]]
    keys <- group$variant_keys[[1]]
    relevant <- controls[
        controls$unit == candidate_unit & variant_key %in% keys
    ]
    relevant <- relevant[
        relevant$model %in% compatible_control_models(candidate_model)
    ]
    control_pids <- sort(unique(relevant$pid))

    list(
        validation = if (length(control_pids) > 0L) "control_match" else "ok",
        control_count = length(control_pids),
        control_pids = collapse_ids(control_pids)
    )
}


validate_candidates <- function(candidates, controls) {
    if (nrow(candidates) == 0L) {
        detail <- copy(candidates)
        detail[, `:=`(
            validation = character(),
            control_count = integer(),
            control_pids = character()
        )]
        summary <- data.table(
            unit = character(),
            model = character(),
            candidate_signature = character(),
            affected_count = integer(),
            affected_pids = character(),
            minimum_PSAP = numeric(),
            validation = character(),
            control_count = integer(),
            control_pids = character()
        )
        return(list(detail = detail, summary = summary))
    }

    groups <- build_candidate_groups(candidates)
    status <- vector("list", nrow(groups))

    for (i in seq_len(nrow(groups))) {
        result <- validate_one_candidate(groups[i], controls)
        status[[i]] <- data.table(
            candidate_id = groups$candidate_id[[i]],
            validation = result$validation,
            control_count = as.integer(result$control_count),
            control_pids = result$control_pids
        )
    }

    status <- rbindlist(status, use.names = TRUE)
    detail <- status[candidates, on = "candidate_id"]
    setcolorder(
        detail,
        c(
            "candidate_id",
            "pid",
            "phenotype",
            "sex",
            "unit",
            "model",
            "variant_key",
            "PSAP",
            "validation",
            "control_count",
            "control_pids"
        )
    )

    candidate_level <- status[groups, on = "candidate_id"]
    summary <- candidate_level[
        ,
        .(
            affected_count = uniqueN(pid),
            affected_pids = collapse_ids(pid),
            minimum_PSAP = minimum_or_na(PSAP),
            validation = validation_priority(validation),
            control_count = uniqueN(unlist(
                strsplit(
                    control_pids[!is.na(control_pids)],
                    ",",
                    fixed = TRUE
                ),
                use.names = FALSE
            )),
            control_pids = combine_id_lists(control_pids)
        ),
        by = .(unit, model, candidate_signature)
    ]

    setorderv(
        detail,
        c("PSAP", "unit", "model", "candidate_id", "variant_key"),
        c(1L, 1L, 1L, 1L, 1L),
        na.last = TRUE
    )
    setorderv(
        summary,
        c("minimum_PSAP", "unit", "model", "candidate_signature"),
        c(1L, 1L, 1L, 1L),
        na.last = TRUE
    )

    list(detail = detail, summary = summary)
}


main <- function() {
    config <- parse_cli(commandArgs(trailingOnly = TRUE))
    setDTthreads(config$threads)

    candidates <- read_candidate_tables(config$candidates)
    controls <- read_controls(config$controls)
    report <- validate_candidates(candidates, controls)

    fwrite(
        report$detail,
        file = config$output_candidates,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )
    fwrite(
        report$summary,
        file = config$output_summary,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        nThread = config$threads
    )

    message("Affected candidate rows: ", nrow(report$detail))
    message("Cohort summary rows: ", nrow(report$summary))
}


if (sys.nframe() == 0) {
    main()
}
