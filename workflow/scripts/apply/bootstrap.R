script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_args) != 1L) {
    stop("Could not determine the executing R script path", call. = FALSE)
}

PSAP_SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_args)))
PSAP_R_LIB_DIR <- normalizePath(
    file.path(PSAP_SCRIPT_DIR, "..", "..", "lib", "r"),
    mustWork = TRUE
)

source(file.path(PSAP_R_LIB_DIR, "common.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "table_contracts.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "variant_identity.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "disease_models.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "genotypes.R"), local = FALSE)
source(file.path(PSAP_R_LIB_DIR, "lookup_contracts.R"), local = FALSE)
