# Candidate/control joins and cohort-level report aggregation.

read_candidate_tables <- function(paths) {
  required <- c("candidate_id", "pid", "phenotype", "sex", "unit", "partition_id", "model", "variant_key", "PSAP")
  tables <- lapply(paths, function(path) {
    candidates <- read_tsv(path)
    require_columns(candidates, required, paste("Candidate table", path))
    candidates[, required, with = FALSE]
  })
  candidates <- rbindlist(tables, use.names = TRUE, fill = FALSE)
  candidates[, PSAP := suppressWarnings(as.numeric(PSAP))]
  if (any(candidates$phenotype != "affected")) abort("Candidate tables contain rows not labelled affected")
  candidates
}


read_controls <- function(path) {
  required <- c("pid", "phenotype", "sex", "unit", "partition_id", "model", "variant_key", "PSAP")
  controls <- read_tsv(path)
  require_columns(controls, required, "Control reference")
  controls <- controls[, required, with = FALSE]
  # ``fread`` infers logical columns for a header-only control table.  Keep
  # join-key types stable when the cohort has no non-reference controls.
  character_columns <- setdiff(required, "PSAP")
  controls[, (character_columns) := lapply(.SD, as.character), .SDcols = character_columns]
  controls[, PSAP := suppressWarnings(as.numeric(PSAP))]
  if (any(controls$phenotype != "unaffected")) abort("Control reference contains rows not labelled unaffected")
  unique(controls, by = c("pid", "unit", "partition_id", "model", "variant_key"))
}


collapse_ids <- function(value) {
  value <- unique(value[!is.na(value) & nzchar(value)])
  if (length(value) == 0L) NA_character_ else paste(sort(value), collapse = ",")
}


combine_id_lists <- function(value) {
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) == 0L) NA_character_ else {
    collapse_ids(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  }
}


validation_priority <- function(value) {
  for (status in unique(value)) validate_status(status, "triage")
  if ("control_match" %in% value) "control_match" else if (all(value == "no_controls")) "no_controls" else "ok"
}


validate_candidate_rows <- function(candidates) {
  inconsistent <- candidates[, .(
    pid_count = uniqueN(pid), unit_count = uniqueN(unit),
    partition_count = uniqueN(partition_id), model_count = uniqueN(model)
  ), by = candidate_id][
    pid_count != 1L | unit_count != 1L | partition_count != 1L | model_count != 1L
  ]
  if (nrow(inconsistent) > 0L) abort("candidate_id values span multiple patients, units, partitions, or models")
  inconsistent_variant <- candidates[, .(variant_count = uniqueN(variant_key)), by = candidate_id][variant_count != 1L]
  if (nrow(inconsistent_variant) > 0L) abort("candidate_id values must identify one variant")
  invisible(candidates)
}


validate_candidates <- function(candidates, controls) {
  if (nrow(candidates) == 0L) {
    detail <- copy(candidates)
    detail[, `:=`(validation = character(), control_count = integer(), control_pids = character())]
    return(list(detail = detail))
  }
  validate_candidate_rows(candidates)
  candidate_keys <- unique(candidates[, .(candidate_id, unit, partition_id, model, variant_key)])
  candidate_controls <- compatible_control_model_table()[
    candidate_keys, on = "model", allow.cartesian = TRUE
  ]
  hits <- controls[candidate_controls,
    on = .(unit, partition_id, variant_key, model = control_model),
    nomatch = NULL,
    allow.cartesian = TRUE,
    .(candidate_id = i.candidate_id, pid = x.pid)
  ]
  hit_summary <- hits[, .(
    control_count = as.integer(uniqueN(pid)),
    control_pids = collapse_ids(pid)
  ), by = candidate_id]
  status <- data.table(
    candidate_id = candidate_keys$candidate_id,
    control_count = NA_integer_,
    control_pids = NA_character_
  )
  status[hit_summary, on = "candidate_id", `:=`(
    control_count = i.control_count,
    control_pids = i.control_pids
  )]
  if (nrow(controls) == 0L) {
    status[, validation := "no_controls"]
  } else {
    status[, validation := fifelse(is.na(control_count), "ok", "control_match")]
  }
  status[is.na(control_count), control_count := 0L]
  setcolorder(status, c("candidate_id", "validation", "control_count", "control_pids"))
  detail <- status[candidates, on = "candidate_id"]
  setcolorder(detail, c(
    "candidate_id", "pid", "phenotype", "sex", "unit", "partition_id", "model", "variant_key",
    "PSAP", "validation", "control_count", "control_pids"
  ))
  setorderv(detail, c("PSAP", "unit", "model", "candidate_id", "variant_key"), c(1L, 1L, 1L, 1L, 1L), na.last = TRUE)
  list(detail = detail)
}
