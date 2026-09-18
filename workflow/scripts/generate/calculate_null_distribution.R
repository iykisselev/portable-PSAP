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
# The only wired absence policy is the explicitly named zero adapter. Invalid
# and excluded panel alleles never reach this input table as observed AF.
scores[is.na(af), af := 0]
probabilities <- observed_af_events(scores$af)
score_max <- as.numeric(snakemake@params[["score_max"]])
bin_count <- as.integer(snakemake@params[["bin_count"]])

calculate_role <- function(role, q) {
  working <- copy(scores)
  working[, q := q]
  result <- lapply(sort(unique(eligible)), function(unit) {
    distribution <- conditional_maximum_distribution(
      working[testing_unit_id == unit, .(score, q)], score_max, bin_count
    )
    list(
      lookup = data.table(
        testing_unit_id = unit,
        role = role,
        bin = seq_len(bin_count),
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
