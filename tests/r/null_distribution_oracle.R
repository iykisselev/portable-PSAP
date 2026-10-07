# Test-only independent oracle for null-distribution event assembly and binning.
#
# Production generation uses the Python implementation; this file is sourced
# only by the R regression suite to provide an independent implementation.

generated_bin_index <- function(score, score_max, bin_width) {
  fixed_width_bin_index(score, score_max, bin_width)
}


observed_af_events <- function(af) {
  af <- as.numeric(af)
  if (any(!is.finite(af) | af < 0 | af > 1)) {
    abort("Observed AF must be finite and within [0, 1]")
  }
  list(het = 2 * af * (1 - af), hom = af^2)
}


zero_unobserved_events <- function(count = 1L) {
  # This is an explicit adapter, not a default-policy decision.  It assigns no
  # population-frequency event probability to panel-absent eligible alleles.
  list(het = rep(0, count), hom = rep(0, count))
}


assemble_calibration_events <- function(
    allele_states,
    allele_membership,
    allele_scores,
    effective_mask,
    unobserved_policy,
    mutation_rates = NULL
) {
  required_columns <- list(
    allele_states = c("CHROM", "POS", "REF", "ALT", "allele_state", "af", "state_reason"),
    allele_membership = c("testing_unit_id", "CHROM", "POS", "REF", "ALT"),
    allele_scores = c("CHROM", "POS", "REF", "ALT", "score"),
    effective_mask = c("CHROM", "START", "END")
  )
  artifacts <- list(
    allele_states = allele_states,
    allele_membership = allele_membership,
    allele_scores = allele_scores,
    effective_mask = effective_mask
  )
  for (artifact_name in names(artifacts)) {
    missing <- setdiff(required_columns[[artifact_name]], names(artifacts[[artifact_name]]))
    if (length(missing)) {
      abort(artifact_name, " is missing required columns: ", paste(missing, collapse = ", "))
    }
  }

  key_columns <- c("CHROM", "POS", "REF", "ALT")
  states <- copy(allele_states)
  memberships <- copy(allele_membership)
  scores <- copy(allele_scores)
  mask <- copy(effective_mask)
  for (table in list(states, memberships, scores)) {
    table[, CHROM := canonical_chromosome(CHROM)]
    table[, POS := as.integer(POS)]
    table[, REF := toupper(as.character(REF))]
    table[, ALT := toupper(as.character(ALT))]
    table[, allele_key := canonical_allele_key(CHROM, POS, REF, ALT)]
  }

  if (anyDuplicated(states$allele_key)) abort("allele_states contains duplicate allele keys")
  if (anyDuplicated(scores$allele_key)) abort("allele_scores contains duplicate allele keys")
  if (anyDuplicated(memberships[, paste(testing_unit_id, allele_key, sep = "\t")])) {
    abort("allele_unit_membership contains duplicate unit/allele pairs")
  }
  if (any(!states$allele_state %chin% c("observed", "absent", "rejected"))) {
    abort("allele_states contains an unsupported state")
  }
  states[, af := suppressWarnings(as.numeric(af))]
  observed <- states$allele_state == "observed"
  if (any(!is.finite(states$af[observed]) | states$af[observed] < 0 | states$af[observed] > 1)) {
    abort("Observed allele states require finite AF within [0, 1]")
  }
  if (any(!is.na(states$af[!observed]))) {
    abort("Absent and rejected allele states must not carry AF")
  }
  rejected <- states$allele_state == "rejected"
  if (any(is.na(states$state_reason[rejected]) | !nzchar(states$state_reason[rejected]))) {
    abort("Rejected allele states require an explicit state_reason")
  }
  scores[, score := suppressWarnings(as.numeric(score))]
  if (any(!is.finite(scores$score))) abort("allele_scores contains a non-finite score")

  mask[, CHROM := canonical_chromosome(CHROM)]
  mask[, `:=`(START = as.integer(START), END = as.integer(END))]
  if (any(is.na(mask$START) | is.na(mask$END) | mask$START < 0 | mask$END <= mask$START)) {
    abort("effective_calibration_mask contains invalid intervals")
  }
  points <- unique(memberships[, .(
    CHROM,
    start = POS,
    end = POS,
    allele_key
  )])
  mask_ranges <- mask[, .(CHROM, start = START + 1L, end = END)]
  setkey(mask_ranges, CHROM, start, end)
  mask_hits <- foverlaps(points, mask_ranges, type = "within", nomatch = 0L)
  if (uniqueN(mask_hits$allele_key) != nrow(points)) {
    abort("allele_unit_membership contains allele(s) outside effective_calibration_mask")
  }

  events <- merge(
    memberships,
    scores[, c(key_columns, "allele_key", "score"), with = FALSE],
    by = c(key_columns, "allele_key"),
    all = FALSE,
    sort = FALSE
  )
  events <- merge(
    events,
    states[, .(allele_key, allele_state, af, state_reason)],
    by = "allele_key",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(events$allele_state))) {
    abort("A scored allele/unit membership has no canonical allele state")
  }

  # Rejected alleles are removed before either absence adapter is allowed to run.
  events <- events[allele_state != "rejected"]
  events[, `:=`(effective_af = af, assignment_source = "observed_af")]
  absent <- events$allele_state == "absent"
  if (identical(unobserved_policy, "zero")) {
    events[absent, `:=`(effective_af = 0, assignment_source = "absent_zero")]
  } else if (identical(unobserved_policy, "mutation_rate")) {
    if (is.null(mutation_rates)) abort("mutation_rate policy requires allele-specific rates")
    rates <- copy(mutation_rates)
    if (!all(c(key_columns, "mu") %chin% names(rates))) {
      abort("Mutation-rate table needs CHROM, POS, REF, ALT, and mu columns")
    }
    rates[, CHROM := canonical_chromosome(CHROM)]
    rates[, POS := as.integer(POS)]
    rates[, REF := toupper(as.character(REF))]
    rates[, ALT := toupper(as.character(ALT))]
    rates[, allele_key := canonical_allele_key(CHROM, POS, REF, ALT)]
    rates[, mu := suppressWarnings(as.numeric(mu))]
    if (anyDuplicated(rates$allele_key) || any(!is.finite(rates$mu) | rates$mu < 0 | rates$mu > 1)) {
      abort("Mutation-rate table has invalid or duplicate allele-specific mu values")
    }
    events <- merge(events, rates[, .(allele_key, mu)], by = "allele_key", all.x = TRUE, sort = FALSE)
    if (any(events$allele_state == "absent" & is.na(events$mu))) {
      abort("Mutation-rate table is missing an assigned mu for panel-absent scored allele(s)")
    }
    events[allele_state == "absent", `:=`(
      effective_af = mu,
      assignment_source = "absent_mutation_rate"
    )]
  } else {
    abort("Unsupported unobserved policy: ", unobserved_policy)
  }
  if (any(!is.finite(events$effective_af))) {
    abort("Event assembly produced an invalid effective allele frequency")
  }
  events[]
}


conditional_maximum_distribution <- function(events, score_max, bin_width) {
  required <- c("score", "q")
  if (!all(required %chin% names(events))) {
    abort("Events need score and q columns")
  }
  score <- as.numeric(events$score)
  q <- as.numeric(events$q)
  if (any(!is.finite(score)) || any(!is.finite(q) | q < 0 | q > 1)) {
    abort("Event scores must be finite and q values must be finite in [0, 1]")
  }
  bin_count <- fixed_width_bin_count(score_max, bin_width)
  bins <- generated_bin_index(score, score_max, bin_width)
  active <- q > 0
  if (!any(active)) {
    return(list(
      status = "zero_event_mass",
      event_mass = 0,
      survival = rep(NA_real_, bin_count),
      max_mass = rep(0, bin_count)
    ))
  }
  q <- q[active]
  bins <- bins[active]
  log_none <- function(selected) {
    values <- q[selected]
    if (any(values == 1)) return(-Inf)
    sum(log1p(-values))
  }
  total_log_none <- log_none(rep(TRUE, length(q)))
  event_mass <- -expm1(total_log_none)
  if (!is.finite(event_mass) || event_mass <= 0) {
    return(list(
      status = "numerical_failure",
      event_mass = event_mass,
      survival = rep(NA_real_, bin_count),
      max_mass = rep(NA_real_, bin_count)
    ))
  }
  survival <- vapply(seq_len(bin_count), function(bin) {
    -expm1(log_none(bins >= bin)) / event_mass
  }, numeric(1))
  max_mass <- vapply(seq_len(bin_count), function(bin) {
    none_above <- exp(log_none(bins > bin))
    none_at_or_above <- exp(log_none(bins >= bin))
    none_above - none_at_or_above
  }, numeric(1))
  list(status = "ok", event_mass = event_mass, survival = survival, max_mass = max_mass)
}
