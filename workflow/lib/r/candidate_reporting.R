# Affected-sample candidate table normalization and deterministic identifiers.

read_candidate_psap <- function(path) {
  psap <- read_tsv(path)
  require_columns(psap, c("Key", "unit", "partition_id", "Dz.Model", "PSAP"), "Affected PSAP table")
  psap[, unit := clean_identifier(unit)]
  psap[, model := clean_identifier(`Dz.Model`)]
  psap[, variant_key := clean_identifier(Key)]
  psap[, partition_id := clean_identifier(partition_id)]
  psap[, PSAP := suppressWarnings(as.numeric(PSAP))]
  unsupported <- setdiff(unique(psap$model[!is.na(psap$model)]), disease_model_names())
  if (length(unsupported) > 0L) abort("Unsupported Dz.Model value(s): %s", paste(unsupported, collapse = ", "))
  invalid <- is.na(psap$unit) | is.na(psap$model) | is.na(psap$variant_key) | is.na(psap$partition_id)
  if (any(invalid)) abort(
    "Affected PSAP table contains %d row(s) with missing unit, model, or variant key",
    sum(invalid)
  )
  unique(psap[, .(unit, partition_id, model, variant_key, PSAP)], by = c("unit", "partition_id", "model", "variant_key", "PSAP"))
}


empty_candidates <- function() {
  data.table(
    candidate_id = character(), pid = character(), phenotype = character(),
    sex = numeric(), unit = character(), partition_id = character(), model = character(),
    variant_key = character(), PSAP = numeric()
  )
}


build_candidates <- function(psap, pid, sex, candidate_id_separator = ":", candidate_id_ordinal_width = 3L) {
  if (!nzchar(candidate_id_separator)) abort("Candidate ID separator must be non-empty")
  if (is.na(candidate_id_ordinal_width) || candidate_id_ordinal_width < 1L) {
    abort("Candidate ID ordinal width must be a positive integer")
  }
  eligible <- psap[model %in% disease_model_names()]
  if (nrow(eligible) == 0L) return(empty_candidates())
  setorderv(eligible, c("model", "unit", "partition_id", "variant_key"))
  eligible[, candidate_number := seq_len(.N), by = .(model, partition_id)]
  eligible[, candidate_id := paste(
    pid, disease_model_tag(model),
    sprintf(paste0("%0", candidate_id_ordinal_width, "d"), candidate_number),
    sep = candidate_id_separator
  )]
  candidates <- eligible[, .(
    candidate_id, pid = as.character(pid), phenotype = "affected",
    sex = sex, unit, partition_id, model, variant_key, PSAP
  )]
  setorderv(candidates, c("PSAP", "candidate_id", "variant_key"), c(1L, 1L, 1L), na.last = TRUE)
  candidates
}


run_candidate_reporting <- function(config) {
  setDTthreads(config$threads)
  candidates <- build_candidates(
    read_candidate_psap(config$input), config$pid, config$sex_value,
    config$candidate_id_separator, config$candidate_id_ordinal_width
  )
  write_tsv(candidates, config$output, threads = config$threads)
  message("Candidate rows written: ", nrow(candidates))
  message("Wrote affected candidates: ", config$output)
}
