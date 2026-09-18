#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

source(file.path(snakemake@scriptdir, "..", "bootstrap.R"))
source(file.path(PSAP_R_LIB_DIR, "null_distribution.R"))

write_distribution <- function(path, rows) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(rows, path, sep = "\t")
}

observed <- fread(snakemake@input[["observed_af"]], showProgress = FALSE)
scores <- fread(snakemake@input[["score_universe"]], showProgress = FALSE)
units <- fread(snakemake@input[["units"]], showProgress = FALSE)
eligible <- units[status == "eligible", testing_unit_id]

observed[, allele_key := paste(CHROM, POS, REF, ALT, sep = ":")]
scores[, allele_key := paste(CHROM, POS, REF, ALT, sep = ":")]
scores <- merge(
  scores[testing_unit_id %chin% eligible],
  observed[, .(allele_key, af)],
  by = "allele_key",
  all.x = TRUE,
  sort = FALSE
)
policy <- as.character(snakemake@params[["unobserved_policy"]])
if (identical(policy, "mutation_rate")) {
  rate_paths <- snakemake@input[["mutation_rate"]]
  if (length(rate_paths) != 1L) abort("mutation_rate policy requires one declared mutation-rate table")
  rates <- fread(rate_paths[[1L]], header = FALSE, showProgress = FALSE)
  if (ncol(rates) < 5L) abort("Mutation-rate table needs CHROM, POS, REF, ALT, and mu columns")
  setnames(rates, names(rates)[1:5], c("CHROM", "POS", "REF", "ALT", "mu"))
  rates[, CHROM := sub("^CHR", "", toupper(trimws(as.character(CHROM))))]
  rates[, allele_key := paste(CHROM, POS, REF, ALT, sep = ":")]
  rates[, mu := as.numeric(mu)]
  if (any(!is.finite(rates$mu) | rates$mu < 0 | rates$mu > 1) || anyDuplicated(rates$allele_key)) {
    abort("Mutation-rate table has invalid or duplicate allele-specific mu values")
  }
  scores <- merge(scores, rates[, .(allele_key, mu)], by = "allele_key", all.x = TRUE, sort = FALSE)
  if (any(is.na(scores$af) & is.na(scores$mu))) {
    abort("Mutation-rate table is missing an assigned mu for panel-absent scored allele(s)")
  }
  scores[is.na(af), af := mu]
} else if (identical(policy, "zero")) {
  # Invalid and excluded panel alleles never reach this input table as observed AF.
  scores[is.na(af), af := 0]
} else {
  abort("Unsupported unobserved policy: ", policy)
}
probabilities <- observed_af_events(scores$af)
score_max <- as.numeric(snakemake@params[["score_max"]])
bin_width <- as.numeric(snakemake@params[["bin_width"]])

calculate_role <- function(role, q) {
  working <- copy(scores)
  working[, q := q]
  result <- lapply(sort(unique(eligible)), function(unit) {
    distribution <- conditional_maximum_distribution(
      working[testing_unit_id == unit, .(score, q)], score_max, bin_width
    )
    list(
      lookup = data.table(
        testing_unit_id = unit,
        role = role,
        bin = seq_len(ceiling(score_max / bin_width)),
        survival = distribution$survival
      ),
      diagnostic = data.table(
        testing_unit_id = unit,
        role = role,
        event_mass = distribution$event_mass,
        status = distribution$status,
        event_count = nrow(working[testing_unit_id == unit])
      )
    )
  })
  list(
    lookup = rbindlist(lapply(result, `[[`, "lookup")),
    diagnostics = rbindlist(lapply(result, `[[`, "diagnostic"))
  )
}

het <- calculate_role("het", probabilities$het)
hom <- calculate_role("hom", probabilities$hom)
write_distribution(snakemake@output[["lookup_het"]], het$lookup)
write_distribution(snakemake@output[["lookup_hom"]], hom$lookup)
write_distribution(snakemake@output[["diagnostics"]], rbindlist(list(het$diagnostics, hom$diagnostics)))
