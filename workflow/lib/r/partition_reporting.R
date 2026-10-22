# Partition-level PSAP selection and per-person/model ACAT aggregation.

json_array <- function(values) {
  jsonlite::toJSON(as.list(values), auto_unbox = TRUE, digits = NA, na = "null")
}


normalize_json_cell <- function(value) {
  value <- trimws(as.character(value))
  quoted <- !is.na(value) & nchar(value) >= 2L & substr(value, 1L, 1L) == '"' & substr(value, nchar(value), nchar(value)) == '"'
  value[quoted] <- gsub('""', '"', substr(value[quoted], 2L, nchar(value[quoted]) - 1L), fixed = TRUE)
  value
}

partition_best_psaps <- function(scored) {
  required <- c("pid", "unit", "model", "partition_id", "PSAP")
  require_columns(scored, required, "Partition scoring results")
  results <- data.table::copy(data.table::as.data.table(scored))
  results[, PSAP := suppressWarnings(as.numeric(PSAP))]
  if (any(!is.finite(results$PSAP) | results$PSAP < 0 | results$PSAP > 1)) {
    abort("Partition scoring results contain an invalid PSAP")
  }
  results[, .(
    PSAP = min(PSAP),
    status = "scored",
    supporting_variant_count = .N
  ), by = .(pid, unit, model, partition_id)]
}


omnibus_partition_psaps <- function(partition_results, configured_partition_ids) {
  required <- c("pid", "unit", "model", "partition_id", "PSAP", "status")
  require_columns(partition_results, required, "Partition summary results")
  if (!is.character(configured_partition_ids) || !length(configured_partition_ids) ||
      anyNA(configured_partition_ids) || anyDuplicated(configured_partition_ids)) {
    abort("Configured partition IDs must be a nonempty unique character vector")
  }
  results <- data.table::copy(data.table::as.data.table(partition_results))
  if (any(!results$partition_id %in% configured_partition_ids)) {
    abort("Partition results contain an undeclared partition ID")
  }
  groups <- unique(results[, .(pid, unit, model)])
  if (!nrow(groups)) return(data.table::data.table(
    pid = character(), unit = character(), model = character(), PSAP = numeric(),
    aggregation_method = character(), aggregation_status = character(),
    configured_partition_count = integer(), contributing_partition_count = integer(),
    contributing_partitions_json = character()
  ))
  groups <- results[groups, on = .(pid, unit, model), {
      available <- x.status == "scored"
      combined <- combine_partition_psaps(x.PSAP[available])
      list(
        combined$PSAP,
        combined$aggregation_status,
        sum(available),
        json_array(x.partition_id[available])
      )
    }, by = .EACHI]
  data.table::setnames(groups, c("V1", "V2", "V3", "V4"),
    c("PSAP", "aggregation_status", "contributing_partition_count", "contributing_partitions_json"))
  groups[, `:=`(
    aggregation_method = "ACAT",
    configured_partition_count = length(configured_partition_ids)
  )]
  data.table::setcolorder(groups, c("pid", "unit", "model", "PSAP", "aggregation_method",
    "aggregation_status", "configured_partition_count", "contributing_partition_count",
    "contributing_partitions_json"))
  groups
}


empty_partition_report <- function() {
  data.table::data.table(
    pid = character(), phenotype = character(), sex = numeric(), testing_unit_id = character(),
    partition_id = character(), model = character(), PSAP = numeric(), selected_score_bin = integer(),
    PSAP_next_bin = numeric(), mid_PSAP_candidate = numeric(), status = character(),
    failed_roles_json = character(), compatible_variant_count = integer(), top_variant_count = integer(),
    top_variant_keys_json = character(), validation = character(), control_count = integer(),
    control_pids = character(), affected_count = integer(), affected_pids = character(), scope_status = character()
  )
}


empty_main_report <- function() {
  data.table::data.table(
    pid = character(), phenotype = character(), sex = numeric(), testing_unit_id = character(),
    model = character(), PSAP = numeric(), aggregation_method = character(), aggregation_status = character(),
    aggregation_input_policy = character(), configured_partition_count = integer(),
    eligible_partition_count = integer(), contributing_partition_count = integer(),
    contributing_partition_ids_json = character(), partition_psaps_json = character(),
    aggregation_inputs_json = character(), partition_statuses_json = character(), raw_zero_count = integer(),
    raw_one_count = integer(), raw_near_one_counts_json = character(), validation = character(),
    control_count = integer(), control_pids = character(), affected_count = integer(), affected_pids = character()
  )
}


json_control_ids <- function(value) {
  value <- value[!is.na(value) & nzchar(value)]
  if (!length(value)) return(json_array(character()))
  json_array(sort(unique(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))))
}


read_partition_statuses <- function(path) {
  statuses <- data.table::as.data.table(read_tsv(path))
  require_columns(statuses, c("testing_unit_id", "partition_id", "status", "failed_roles_json"), "Certified unit/partition statuses")
  statuses[, `:=`(testing_unit_id = clean_identifier(testing_unit_id), partition_id = clean_identifier(partition_id), failed_roles_json = normalize_json_cell(failed_roles_json))]
  if (anyNA(statuses$testing_unit_id) | anyNA(statuses$partition_id) || anyDuplicated(statuses, by = c("testing_unit_id", "partition_id"))) {
    abort("Certified unit/partition statuses must have unique nonmissing unit/partition keys")
  }
  statuses
}


compatible_evidence_universe <- function(paths, pids, sexes) {
  if (length(paths) != length(pids) || length(paths) != length(sexes)) {
    abort("Affected evidence paths, PIDs, and sexes must have equal lengths")
  }
  evidence <- data.table::rbindlist(lapply(seq_along(paths), function(index) {
    rows <- data.table::as.data.table(read_tsv(paths[[index]]))
    require_columns(rows, c("unit", "partition_id", "CHROM", "gt_class", "partition_status"), "Affected membership evidence")
    rows[, `:=`(pid = as.character(pids[[index]]), phenotype = "affected", sex = as.numeric(sexes[[index]]),
      unit = clean_identifier(unit), partition_id = clean_identifier(partition_id),
      gt_class = tolower(trimws(as.character(gt_class))))]
    rows
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(evidence)) return(data.table::data.table(pid = character(), phenotype = character(), sex = numeric(), unit = character(), model = character()))
  evidence[, scope := ifelse(canonical_chromosome(CHROM) == "X", "non_PAR_chrX", "autosomal")]
  evidence <- evidence[scope %in% c("autosomal", "non_PAR_chrX")]
  definitions <- data.table::as.data.table(DISEASE_MODELS)[, .(model, gt_class = genotype_class, chromosome_scope)]
  result <- merge(evidence, definitions, by = "gt_class", allow.cartesian = TRUE)
  result <- result[scope == chromosome_scope]
  unique(result[, .(pid, phenotype, sex, unit, model)], by = c("pid", "unit", "model"))
}


partitioned_cohort_reports <- function(candidates, configured_partition_ids,
                                       available_evidence = NULL, affected_pids = NULL,
                                       affected_sexes = NULL, partition_statuses = NULL) {
  required <- c("pid", "phenotype", "sex", "unit", "partition_id", "model", "variant_key", "PSAP", "selected_score_bin", "next_PSAP", "validation", "control_count", "control_pids")
  require_columns(candidates, required, "Partitioned affected candidates")
  scored <- data.table::copy(data.table::as.data.table(candidates))
  scored[, `:=`(PSAP = as.numeric(PSAP), selected_score_bin = as.integer(selected_score_bin), next_PSAP = as.numeric(next_PSAP))]
  if (any(!is.finite(scored$PSAP) | scored$PSAP < 0 | scored$PSAP > 1)) abort("Partitioned affected candidates contain an invalid PSAP")
  if (is.null(available_evidence)) {
    universe <- unique(scored[, .(pid, phenotype, sex, unit, model)])
  } else {
    universe <- compatible_evidence_universe(available_evidence, affected_pids, affected_sexes)
  }
  if (is.null(partition_statuses)) {
    statuses <- unique(scored[, .(testing_unit_id = unit, partition_id, status = "eligible", failed_roles_json = "[]")])
  } else {
    statuses <- read_partition_statuses(partition_statuses)
  }
  if (!nrow(universe)) return(list(main = empty_main_report(), partitions = replicate(length(configured_partition_ids), empty_partition_report(), simplify = FALSE), detail = empty_partition_report()))

  rows <- vector("list", nrow(universe) * length(configured_partition_ids))
  cursor <- 0L
  for (group_index in seq_len(nrow(universe))) {
    group <- universe[group_index]
    for (partition_id in configured_partition_ids) {
      cursor <- cursor + 1L
      target_partition_id <- partition_id
      scope <- statuses[testing_unit_id == group$unit & partition_id == target_partition_id]
      scope_status <- if (!nrow(scope)) "unit_absent" else scope$status[[1L]]
      failed_roles <- if (!nrow(scope)) "[]" else scope$failed_roles_json[[1L]]
      hit <- scored[pid == group$pid & unit == group$unit & model == group$model & partition_id == target_partition_id]
      if (nrow(hit)) {
        best_psap <- min(hit$PSAP)
        top <- hit[PSAP == best_psap]
        selected_bin <- max(top$selected_score_bin)
        next_psap <- top[selected_score_bin == selected_bin, next_PSAP[[1L]]]
        status <- "scored"
        failed_roles <- "[]"
        validation <- validation_priority(top$validation)
        control_values <- top$control_pids
        control_count <- length(unique(unlist(strsplit(control_values[!is.na(control_values) & nzchar(control_values)], ",", fixed = TRUE), use.names = FALSE)))
        top_keys <- json_array(sort(unique(top$variant_key)))
        compatible_count <- data.table::uniqueN(hit$variant_key)
        top_count <- data.table::uniqueN(top$variant_key)
      } else {
        best_psap <- NA_real_; selected_bin <- NA_integer_; next_psap <- NA_real_
        status <- if (scope_status == "eligible") "no_compatible_variant" else scope_status
        validation <- "not_evaluated"; control_values <- character(); control_count <- 0L
        top_keys <- json_array(character()); compatible_count <- 0L; top_count <- 0L
      }
      rows[[cursor]] <- data.table::data.table(
        pid = group$pid, phenotype = group$phenotype, sex = group$sex, testing_unit_id = group$unit,
        partition_id = partition_id, model = group$model, PSAP = best_psap,
        selected_score_bin = selected_bin, PSAP_next_bin = next_psap,
        mid_PSAP_candidate = if (is.na(best_psap)) NA_real_ else (best_psap + next_psap) / 2,
        status = status, failed_roles_json = failed_roles,
        compatible_variant_count = as.integer(compatible_count), top_variant_count = as.integer(top_count),
        top_variant_keys_json = top_keys, validation = validation, control_count = as.integer(control_count),
        control_pids = json_control_ids(control_values), affected_count = 1L,
        affected_pids = json_array(group$pid), scope_status = scope_status
      )
    }
  }
  detail <- data.table::rbindlist(rows, use.names = TRUE)
  main <- detail[, {
    ordered <- .SD[match(configured_partition_ids, partition_id)]
    available <- ordered$status == "scored"
    combined <- combine_partition_psaps(ordered$PSAP[available])
    raw_pairs <- lapply(seq_along(configured_partition_ids), function(i) list(id = configured_partition_ids[[i]], value = ordered$PSAP[[i]]))
    status_pairs <- lapply(seq_along(configured_partition_ids), function(i) list(id = configured_partition_ids[[i]], status = ordered$status[[i]], failed_roles = jsonlite::fromJSON(ordered$failed_roles_json[[i]])))
    controls <- ordered$control_pids[available]
    list(
      PSAP = combined$PSAP, aggregation_method = "ACAT", aggregation_status = combined$aggregation_status,
      aggregation_input_policy = "raw_survival_v1", configured_partition_count = length(configured_partition_ids),
      eligible_partition_count = sum(ordered$scope_status == "eligible"), contributing_partition_count = sum(available),
      contributing_partition_ids_json = json_array(configured_partition_ids[available]),
      partition_psaps_json = jsonlite::toJSON(raw_pairs, auto_unbox = TRUE, digits = NA, na = "null"),
      aggregation_inputs_json = json_array(combined$aggregation_inputs),
      partition_statuses_json = jsonlite::toJSON(status_pairs, auto_unbox = TRUE, digits = NA, na = "null"),
      raw_zero_count = sum(ordered$PSAP[available] == 0), raw_one_count = sum(ordered$PSAP[available] == 1),
      raw_near_one_counts_json = jsonlite::toJSON(list(`0.99` = sum(ordered$PSAP[available] >= 0.99), `0.999999` = sum(ordered$PSAP[available] >= 0.999999), `0.999999999` = sum(ordered$PSAP[available] >= 0.999999999)), auto_unbox = TRUE, digits = NA, na = "null"),
      validation = if (!any(available)) "not_evaluated" else validation_priority(ordered$validation[available]),
      control_count = length(unique(unlist(lapply(controls, function(cell) jsonlite::fromJSON(cell)), use.names = FALSE))),
      control_pids = json_array(sort(unique(unlist(lapply(controls, function(cell) jsonlite::fromJSON(cell)), use.names = FALSE)))),
      affected_count = 1L, affected_pids = json_array(pid[[1L]])
    )
  }, by = .(pid, phenotype, sex, testing_unit_id, model)]
  data.table::setorderv(detail, c("PSAP", "testing_unit_id", "model", "pid", "partition_id"), na.last = TRUE)
  data.table::setorderv(main, c("PSAP", "testing_unit_id", "model", "pid"), na.last = TRUE)
  partitions <- lapply(configured_partition_ids, function(id) detail[partition_id == id])
  list(main = main, partitions = partitions, detail = detail)
}
