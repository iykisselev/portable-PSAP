#!/usr/bin/env Rscript

# Construct patient-specific candidates from one affected sample's PSAP output.

suppressPackageStartupMessages({
    library(data.table)
})

if (exists("snakemake", inherits = FALSE)) {
    source(file.path(snakemake@scriptdir, "bootstrap.R"))
}

options(scipen = 999)


workflow_config <- function() {
    pid <- as.character(snakemake@params[["pid"]])
    if (!nzchar(trimws(pid))) abort("Sample ID must be non-empty")
    list(
        input = snakemake@input[["psap"]],
        pid = pid,
        sex_value = parse_sex(as.character(snakemake@params[["sex"]])),
        candidate_id_separator = as.character(
            snakemake@params[["candidate_id_separator"]]
        ),
        candidate_id_ordinal_width = as.integer(
            snakemake@params[["candidate_id_ordinal_width"]]
        ),
        output = snakemake@output[["candidates"]],
        threads = snakemake@threads
    )
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


build_candidates <- function(
    psap,
    pid,
    sex,
    candidate_id_separator = ":",
    candidate_id_ordinal_width = 3L
) {
    if (!nzchar(candidate_id_separator)) {
        abort("Candidate ID separator must be non-empty")
    }
    if (is.na(candidate_id_ordinal_width) || candidate_id_ordinal_width < 1L) {
        abort("Candidate ID ordinal width must be a positive integer")
    }
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
            sprintf(
                paste0("%0", candidate_id_ordinal_width, "d"),
                candidate_number
            ),
            sep = candidate_id_separator
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


run <- function(config) {
    setDTthreads(config$threads)

    psap <- read_psap(config$input)
    candidates <- build_candidates(
        psap,
        config$pid,
        config$sex_value,
        config$candidate_id_separator,
        config$candidate_id_ordinal_width
    )

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


if (exists("snakemake", inherits = FALSE)) {
    run(workflow_config())
}
