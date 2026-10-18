# Partition-level PSAP selection and per-person/model ACAT aggregation.

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
        jsonlite::toJSON(x.partition_id[available], auto_unbox = TRUE, digits = NA, na = "null")
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
