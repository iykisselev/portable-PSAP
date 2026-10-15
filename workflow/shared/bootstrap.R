# Load the workflow R library when a Snakemake script starts.
PSAP_SCRIPT_DIR <- normalizePath(snakemake@scriptdir, mustWork = TRUE)
psap_find_workflow_root <- function(script_dir) {
  candidate <- script_dir
  repeat {
    if (dir.exists(file.path(candidate, "lib", "r"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop("Could not find workflow/lib/r above script directory: ", script_dir)
}
PSAP_R_LIB_DIR <- file.path(psap_find_workflow_root(PSAP_SCRIPT_DIR), "lib", "r")

# Source the complete workflow library deliberately: each R adapter has the same
# stable namespace, and the shared-library dependency boundary remains explicit here.
source(file.path(PSAP_R_LIB_DIR, "common.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "model_catalog.R"), local = FALSE)
PSAP_MODEL_CATALOG <- read_model_catalog(file.path(dirname(dirname(PSAP_R_LIB_DIR)), "model_catalog.json"))
source(file.path(PSAP_R_LIB_DIR, "acat.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "table_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "genomic_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "variant_identity.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "disease_models.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "genotypes.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "bundle_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "lookup_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "psap_scoring.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "cohort_reporting.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "control_reference.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "candidate_reporting.R"), local = FALSE)
