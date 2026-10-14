# Versioned generated-bundle contract constants shared by R lookup consumers.

CURRENT_GENERATED_BUNDLE_CONTRACT <- "generated_bundle_v2"
INDEL_GENERATED_BUNDLE_CONTRACT <- "generated_bundle_v3"
SUPPORTED_GENERATED_BUNDLE_CONTRACTS <- c(
  "generated_bundle_v1",
  CURRENT_GENERATED_BUNDLE_CONTRACT,
  INDEL_GENERATED_BUNDLE_CONTRACT
)
FIXED_WIDTH_BINNING_ALGORITHM <- "fixed_width_decimal_v1"

is_supported_generated_bundle_contract <- function(contract) {
  as.character(contract) %chin% SUPPORTED_GENERATED_BUNDLE_CONTRACTS
}
