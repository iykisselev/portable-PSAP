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
  if (model == "DOM-het") disease_model_names() else model
}
