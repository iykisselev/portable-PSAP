# Partition-level PSAP selection and per-person/model ACAT aggregation.

json_array <- function(values) {
  jsonlite::toJSON(as.list(values), auto_unbox = TRUE, digits = NA, na = "null")
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


partitioned_cohort_reports <- function(candidates, configured_partition_ids) {
  required <- c("pid", "phenotype", "sex", "unit", "partition_id", "model", "variant_key", "PSAP", "validation", "control_count", "control_pids")
  require_columns(candidates, required, "Partitioned affected candidates")
  scored <- data.table::copy(data.table::as.data.table(candidates))
  scored[, PSAP := as.numeric(PSAP)]
  best <- scored[, .(
    PSAP = min(PSAP),
    supporting_variant_count = .N,
    supporting_variant_keys_json = json_array(sort(unique(variant_key[PSAP == min(PSAP)]))),
    validation = validation_priority(validation),
    control_count = uniqueN(unlist(strsplit(control_pids[!is.na(control_pids)], ",", fixed = TRUE), use.names = FALSE)),
    control_pids = combine_id_lists(control_pids)
  ), by = .(pid, phenotype, sex, unit, partition_id, model)]
  partition_outputs <- lapply(configured_partition_ids, function(partition_id) {
    target_partition_id <- partition_id
    result <- best[partition_id == target_partition_id]
    data.table::setorder(result, PSAP, pid, unit, model)
    result
  })
  main <- best[, {
    ordered <- .SD[match(configured_partition_ids, partition_id)]
    available <- !is.na(ordered$PSAP)
    combined <- combine_partition_psaps(ordered$PSAP[available])
    pairs <- lapply(seq_along(configured_partition_ids), function(i) list(
      id = configured_partition_ids[[i]],
      value = if (available[[i]]) ordered$PSAP[[i]] else NA_real_
    ))
    list(
      PSAP = combined$PSAP,
      aggregation_method = "ACAT",
      aggregation_status = combined$aggregation_status,
      aggregation_input_policy = "raw_survival_v1",
      configured_partition_count = length(configured_partition_ids),
      eligible_partition_count = sum(available),
      contributing_partition_count = sum(available),
      contributing_partition_ids_json = json_array(configured_partition_ids[available]),
      partition_psaps_json = jsonlite::toJSON(pairs, auto_unbox = TRUE, digits = NA, na = "null"),
      aggregation_inputs_json = json_array(combined$aggregation_inputs),
      partition_statuses_json = jsonlite::toJSON(lapply(configured_partition_ids, function(id) list(id = id, status = if (id %in% ordered$partition_id) "scored" else "no_compatible_variant")), auto_unbox = TRUE, digits = NA, na = "null"),
      raw_zero_count = sum(ordered$PSAP[available] == 0),
      raw_one_count = sum(ordered$PSAP[available] == 1),
      raw_near_one_counts_json = jsonlite::toJSON(list(`0.99` = sum(ordered$PSAP[available] >= 0.99)), auto_unbox = TRUE, digits = NA, na = "null"),
      validation = validation_priority(ordered$validation[available]),
      control_count = uniqueN(unlist(strsplit(ordered$control_pids[available & !is.na(ordered$control_pids)], ",", fixed = TRUE), use.names = FALSE)),
      control_pids = combine_id_lists(ordered$control_pids[available]),
      affected_count = 1L, affected_pids = json_array(pid[[1]])
    )
  }, by = .(pid, phenotype, sex, unit, model)]
  data.table::setorderv(main, c("PSAP", "pid", "unit", "model"), na.last = TRUE)
  list(main = main, partitions = partition_outputs, detail = best)
}
