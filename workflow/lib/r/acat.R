# Aggregated Cauchy Association Test helpers for partition-level PSAP values.

validate_acat_inputs <- function(values, label = "ACAT inputs") {
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0 | values > 1)) {
    abort(paste0(label, " must be finite numeric values in [0, 1]"))
  }
  invisible(values)
}


acat_interior <- function(values) {
  # This is deliberately only the continuous, open-interval combiner. Endpoint
  # policy belongs in combine_partition_psaps(), where it remains auditable.
  validate_acat_inputs(values)
  if (!length(values) || any(values <= 0 | values >= 1)) {
    abort("acat_interior requires one or more values strictly between 0 and 1")
  }
  transforms <- numeric(length(values))
  non_half <- values != 0.5
  # tanpi(0.5) is NaN in base R, though its cotangent limit is exactly zero.
  transforms[non_half] <- 1 / tanpi(values[non_half])
  statistic <- sum(transforms) / length(transforms)
  result <- stats::pcauchy(statistic, lower.tail = FALSE)
  if (!is.finite(result) || result < 0 || result > 1) {
    abort("ACAT numerical failure")
  }
  result
}


combine_partition_psaps <- function(raw_psap, aggregation_inputs = raw_psap,
                                    input_policy = c("raw", "mid_p")) {
  # raw_psap is reported as-is in partition summaries. aggregation_inputs is
  # separately supplied for the measured mid-p policy experiment.
  input_policy <- match.arg(input_policy)
  validate_acat_inputs(raw_psap, "Raw partition PSAPs")
  validate_acat_inputs(aggregation_inputs, "ACAT inputs")
  if (length(raw_psap) != length(aggregation_inputs)) {
    abort("Raw partition PSAPs and ACAT inputs must have the same length")
  }
  k <- length(raw_psap)
  if (!k) return(list(PSAP = NA_real_, aggregation_status = "no_contributing_partitions"))

  any_zero <- any(raw_psap == 0)
  any_one <- any(raw_psap == 1)
  if (any_zero) {
    return(list(
      PSAP = 0,
      aggregation_status = if (any_one) "zero_and_one_partition_psap" else "zero_partition_psap"
    ))
  }
  if (k == 1L) {
    return(list(
      PSAP = raw_psap[[1L]],
      aggregation_status = if (raw_psap[[1L]] == 1) "single_partition_one" else "single_partition"
    ))
  }
  if (any(aggregation_inputs == 1)) {
    return(list(PSAP = NA_real_, aggregation_status = "upper_endpoint_unresolved"))
  }
  # A zero aggregation input that did not originate as an exact raw PSAP is
  # invalid under the current policies and must not be silently transformed.
  if (any(aggregation_inputs == 0)) abort("ACAT inputs contain an unexplained zero")
  status <- if (input_policy == "mid_p") {
    if (any_one) "mid_p_exact_one_adjusted" else "ok_mid_p"
  } else {
    "ok"
  }
  list(PSAP = acat_interior(aggregation_inputs), aggregation_status = status)
}
