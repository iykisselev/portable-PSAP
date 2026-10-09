DISEASE_MODELS <- data.frame(
  model = c(
    "DOM-het",
    "REC-hom",
    "X-linked-hom",
    "X-linked-het",
    "X-linked-hem"
  ),
  genotype_class = c(
    "het",
    "hom",
    "hom_female",
    "het_female",
    "hem"
  ),
  lookup = c("het", "hom", "hom", "het", "hem_chrX"),
  tag = c("DOM", "REC_HOM", "XL_HOM", "XL_HET", "XL_HEM"),
  is_hem = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  required_chromosome = c(NA, NA, NA, NA, "X"),
  compatible_models = I(list(
    c("DOM-het", "REC-hom", "X-linked-hom", "X-linked-het", "X-linked-hem"),
    "REC-hom",
    "X-linked-hom",
    "X-linked-het",
    "X-linked-hem"
  )),
  stringsAsFactors = FALSE
)


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
