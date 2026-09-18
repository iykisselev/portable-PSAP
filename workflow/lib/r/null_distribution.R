# Pure event adapters and stable conditional maximum-score distributions.

generated_bin_index <- function(score, score_max, bin_width) {
  if (!is.finite(score_max) || score_max <= 0 || !is.finite(bin_width) || bin_width <= 0) {
    abort("score_max and bin_width must be positive")
  }
  bin_count <- ceiling(score_max / bin_width)
  clamped <- pmin(pmax(as.numeric(score), 0), score_max)
  pmin(floor(clamped / bin_width) + 1L, bin_count)
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
  bin_count <- ceiling(score_max / bin_width)
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
