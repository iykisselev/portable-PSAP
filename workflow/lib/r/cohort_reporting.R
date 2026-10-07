# Candidate/control joins and cohort-level report aggregation.

read_candidate_tables <- function(paths) {
  required <- c("candidate_id", "pid", "phenotype", "sex", "unit", "model", "variant_key", "PSAP")
  tables <- lapply(paths, function(path) {
    candidates <- fread(path, sep = "\t", header = TRUE, quote = "", data.table = TRUE, showProgress = FALSE)
    require_columns(candidates, required, paste("Candidate table", path))
    candidates[, required, with = FALSE]
  })
  candidates <- rbindlist(tables, use.names = TRUE, fill = FALSE)
  candidates[, PSAP := suppressWarnings(as.numeric(PSAP))]
  if (any(candidates$phenotype != "affected")) abort("Candidate tables contain rows not labelled affected")
  candidates
}


read_controls <- function(path) {
  required <- c("pid", "phenotype", "sex", "unit", "model", "variant_key", "PSAP")
  controls <- fread(path, sep = "\t", header = TRUE, quote = "", data.table = TRUE, showProgress = FALSE)
  require_columns(controls, required, "Control reference")
  controls <- controls[, required, with = FALSE]
  # ``fread`` infers logical columns for a header-only control table.  Keep
  # join-key types stable when the cohort has no non-reference controls.
  character_columns <- setdiff(required, "PSAP")
  controls[, (character_columns) := lapply(.SD, as.character), .SDcols = character_columns]
  controls[, PSAP := suppressWarnings(as.numeric(PSAP))]
  if (any(controls$phenotype != "unaffected")) abort("Control reference contains rows not labelled unaffected")
  unique(controls, by = c("pid", "unit", "model", "variant_key"))
}


minimum_or_na <- function(value) {
  if (length(value) == 0L || all(is.na(value))) NA_real_ else min(value, na.rm = TRUE)
}


collapse_ids <- function(value) {
  value <- unique(value[!is.na(value) & nzchar(value)])
  if (length(value) == 0L) NA_character_ else paste(sort(value), collapse = ",")
}


combine_id_lists <- function(value) {
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) == 0L) NA_character_ else collapse_ids(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
}


validation_priority <- function(value) {
  if ("control_match" %in% value) "control_match" else if (all(value == "no_controls")) "no_controls" else "ok"
}


build_candidate_groups <- function(candidates) {
  inconsistent <- candidates[, .(pid_count = uniqueN(pid), unit_count = uniqueN(unit), model_count = uniqueN(model)), by = candidate_id][pid_count != 1L | unit_count != 1L | model_count != 1L]
  if (nrow(inconsistent) > 0L) abort("candidate_id values span multiple patients, units, or models")
  groups <- candidates[, .(pid = pid[[1]], unit = unit[[1]], model = model[[1]], variant_keys = list(sort(unique(variant_key))), PSAP = minimum_or_na(PSAP)), by = candidate_id]
  groups[, variant_count := lengths(variant_keys)]
  if (any(groups$variant_count != 1L)) abort("Candidates must contain exactly one distinct variant")
  groups[, candidate_signature := vapply(variant_keys, paste, collapse = "|", FUN.VALUE = character(1))]
  groups
}


empty_cohort_summary <- function() {
  data.table(unit = character(), model = character(), candidate_signature = character(), affected_count = integer(), affected_pids = character(), minimum_PSAP = numeric(), validation = character(), control_count = integer(), control_pids = character())
}


validate_candidates <- function(candidates, controls) {
  if (nrow(candidates) == 0L) {
    detail <- copy(candidates)
    detail[, `:=`(validation = character(), control_count = integer(), control_pids = character())]
    return(list(detail = detail, summary = empty_cohort_summary()))
  }
  groups <- build_candidate_groups(candidates)
  candidate_keys <- groups[, .(
    candidate_id,
    unit,
    model,
    variant_key = vapply(variant_keys, function(keys) keys[[1L]], character(1))
  )]
  models <- disease_model_names()
  compatible_models <- rbind(
    data.table(model = "DOM-het", control_model = models),
    data.table(
      model = setdiff(models, "DOM-het"),
      control_model = setdiff(models, "DOM-het")
    )
  )
  candidate_controls <- compatible_models[candidate_keys, on = "model", allow.cartesian = TRUE]
  hits <- controls[candidate_controls,
    on = .(unit, variant_key, model = control_model),
    nomatch = NULL,
    allow.cartesian = TRUE,
    .(candidate_id = i.candidate_id, pid = x.pid)
  ]
  hit_summary <- hits[, .(
    control_count = as.integer(uniqueN(pid)),
    control_pids = collapse_ids(pid)
  ), by = candidate_id]
  status <- data.table(
    candidate_id = groups$candidate_id,
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
  setcolorder(detail, c("candidate_id", "pid", "phenotype", "sex", "unit", "model", "variant_key", "PSAP", "validation", "control_count", "control_pids"))
  candidate_level <- status[groups, on = "candidate_id"]
  summary <- candidate_level[, .(
    affected_count = uniqueN(pid), affected_pids = collapse_ids(pid), minimum_PSAP = minimum_or_na(PSAP),
    validation = validation_priority(validation),
    control_count = uniqueN(unlist(strsplit(control_pids[!is.na(control_pids)], ",", fixed = TRUE), use.names = FALSE)),
    control_pids = combine_id_lists(control_pids)
  ), by = .(unit, model, candidate_signature)]
  setorderv(detail, c("PSAP", "unit", "model", "candidate_id", "variant_key"), c(1L, 1L, 1L, 1L, 1L), na.last = TRUE)
  setorderv(summary, c("minimum_PSAP", "unit", "model", "candidate_signature"), c(1L, 1L, 1L, 1L), na.last = TRUE)
  list(detail = detail, summary = summary)
}


run_cohort_reporting <- function(config) {
  setDTthreads(config$threads)
  report <- validate_candidates(read_candidate_tables(config$candidates), read_controls(config$controls))
  fwrite(report$detail, file = config$output_candidates, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, nThread = config$threads)
  fwrite(report$summary, file = config$output_summary, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, nThread = config$threads)
  message("Affected candidate rows: ", nrow(report$detail))
  message("Cohort summary rows: ", nrow(report$summary))
}
