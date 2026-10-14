# Preserve catalog order for model selection, reporting, and control matching.
DISEASE_MODELS <- do.call(rbind, lapply(PSAP_MODEL_CATALOG$models, function(model) {
  data.frame(
    model = model$id,
    genotype_class = model$genotype_class,
    lookup = model$role,
    tag = model$tag,
    is_hem = model$requires_hemizygote,
    required_chromosome = if (model$scope == "non_PAR_chrX") "X" else NA_character_,
    compatible_models = I(list(model$compatible_control_models)),
    chromosome_scope = model$scope,
    stringsAsFactors = FALSE
  )
}))


disease_model_names <- function() {
  DISEASE_MODELS$model
}


disease_model_tag <- function(model) {
  definitions <- DISEASE_MODELS[match(model, DISEASE_MODELS$model), , drop = FALSE]
  definitions$tag
}


compatible_control_models <- function(model) {
  definition <- DISEASE_MODELS[match(model, DISEASE_MODELS$model), , drop = FALSE]
  if (nrow(definition) != 1L) abort("Unsupported disease model: %s", model)
  definition$compatible_models[[1L]]
}


standard_disease_models <- function() {
  DISEASE_MODELS[!DISEASE_MODELS$is_hem, , drop = FALSE]
}


hem_disease_model <- function() {
  definition <- DISEASE_MODELS[DISEASE_MODELS$is_hem, , drop = FALSE]
  if (nrow(definition) != 1L) abort("Disease model contract requires one hemizygous model")
  definition
}


compatible_control_model_table <- function() {
  rbindlist(lapply(seq_len(nrow(DISEASE_MODELS)), function(index) {
    definition <- DISEASE_MODELS[index, , drop = FALSE]
    data.table(
      model = definition$model,
      control_model = compatible_control_models(definition$model)
    )
  }))
}
