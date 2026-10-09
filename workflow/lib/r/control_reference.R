# Reusable extraction and table joins for unaffected control evidence.

empty_reference <- function() {
  data.table(pid = character(), phenotype = character(), sex = numeric(), unit = character(), model = character(), variant_key = character(), PSAP = numeric())
}

make_evidence <- function(variants, row_filter, model) {
  selected <- variants[row_filter]
  if (nrow(selected) == 0L) return(NULL)
  selected[, .(unit, model = model, variant_key)]
}

derive_control_evidence <- function(path, pid, sex, hemizygote_model) {
  variants <- read_tsv(path)
  require_columns(variants, c("Key", "gt_class"), sprintf("Control annotation table for %s", pid))
  variants[, variant_key := clean_identifier(Key)]
  variants[, gt_class := tolower(clean_identifier(gt_class))]
  require_columns(variants, "unit", sprintf("Control annotation table for %s", pid))
  variants[, unit := clean_identifier(unit)]
  invalid <- is.na(variants$variant_key) | is.na(variants$gt_class) | is.na(variants$unit)
  if (any(invalid)) abort("Control %s contains %d row(s) with missing key, genotype class, or unit", pid, sum(invalid))
  definitions <- standard_disease_models()
  evidence <- lapply(seq_len(nrow(definitions)), function(index) {
    definition <- definitions[index, ]
    make_evidence(variants, variants$gt_class == definition$genotype_class, definition$model)
  })
  if (hemizygote_model && any(variants$gt_class == "hem")) {
    require_columns(variants, "CHROM", sprintf("Hemizygous control annotation table for %s", pid))
    chromosome <- canonical_chromosome(variants$CHROM)
    hem_definition <- hem_disease_model()
    evidence <- c(evidence, list(make_evidence(variants, variants$gt_class == hem_definition$genotype_class & chromosome == hem_definition$required_chromosome, hem_definition$model)))
  }
  evidence <- Filter(Negate(is.null), evidence)
  if (length(evidence) == 0L) return(empty_reference())
  output <- unique(rbindlist(evidence, use.names = TRUE), by = c("unit", "model", "variant_key"))
  output[, `:=`(pid = as.character(pid), phenotype = "unaffected", sex = sex, PSAP = NA_real_)]
  setcolorder(output, c("pid", "phenotype", "sex", "unit", "model", "variant_key", "PSAP"))
  output
}


read_psap_scores <- function(path, pid) {
  scores <- read_tsv(path)
  require_columns(scores, c("Key", "unit", "Dz.Model", "PSAP"), sprintf("Control PSAP table for %s", pid))
  scores[, unit := clean_identifier(unit)]
  scores[, model := clean_identifier(`Dz.Model`)]
  scores[, variant_key := clean_identifier(Key)]
  scores[, PSAP := suppressWarnings(as.numeric(PSAP))]
  scores <- scores[!is.na(unit) & !is.na(model) & !is.na(variant_key), .(unit, model, variant_key, PSAP)]
  scores[, .(PSAP = if (all(is.na(PSAP))) NA_real_ else min(PSAP, na.rm = TRUE)), by = .(unit, model, variant_key)]
}


attach_psap <- function(evidence, scores) {
  if (nrow(evidence) == 0L || nrow(scores) == 0L) return(evidence)
  evidence[, PSAP := NULL]
  scores[evidence, on = c("unit", "model", "variant_key")][, .(pid, phenotype, sex, unit, model, variant_key, PSAP)]
}


run_control_reference <- function(config) {
  setDTthreads(config$threads)
  if (length(config$input) == 0L) {
    reference <- empty_reference()
  } else {
    controls <- vector("list", length(config$input))
    for (i in seq_along(config$input)) {
      evidence <- derive_control_evidence(config$input[[i]], config$pid[[i]], config$sex_value[[i]], config$hemizygote_model)
      controls[[i]] <- attach_psap(evidence, read_psap_scores(config$psap_input[[i]], config$pid[[i]]))
    }
    reference <- rbindlist(controls, use.names = TRUE, fill = FALSE)
  }
  if (config$deduplicate && nrow(reference) > 0L) {
    before <- nrow(reference)
    reference <- unique(reference, by = names(reference))
    message("Exact duplicate rows removed: ", before - nrow(reference))
  }
  if (nrow(reference) > 0L) setorderv(reference, c("unit", "model", "variant_key", "pid"), c(1L, 1L, 1L, 1L))
  write_tsv(reference, config$output, threads = config$threads)
  message("Unaffected individuals: ", length(config$input))
  message("Control evidence rows written: ", nrow(reference))
  message("Wrote control reference: ", config$output)
}
