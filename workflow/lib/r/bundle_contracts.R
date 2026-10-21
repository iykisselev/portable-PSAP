# Versioned generated-bundle contract constants shared by R lookup consumers.

PARTITIONED_NULL_BUNDLE_CONTRACT <- "partitioned_null_bundle_v1"
SUPPORTED_GENERATED_BUNDLE_CONTRACTS <- PARTITIONED_NULL_BUNDLE_CONTRACT
FIXED_WIDTH_BINNING_ALGORITHM <- "fixed_width_decimal_v1"

is_supported_generated_bundle_contract <- function(contract) {
  as.character(contract) %chin% SUPPORTED_GENERATED_BUNDLE_CONTRACTS
}
