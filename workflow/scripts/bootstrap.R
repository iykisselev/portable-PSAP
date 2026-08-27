# Locate the workflow R library when invoked by Snakemake
PSAP_SCRIPT_DIR <- normalizePath(snakemake@scriptdir, mustWork = TRUE) 
PSAP_R_LIB_DIR <- normalizePath(
  file.path(PSAP_SCRIPT_DIR, "..", "lib", "r"),
  mustWork = TRUE
)

source(file.path(PSAP_R_LIB_DIR, "common.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "table_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "variant_identity.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "disease_models.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "genotypes.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "lookup_contracts.R"), local = FALSE)
