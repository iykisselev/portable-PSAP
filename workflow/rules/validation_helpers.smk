from pathlib import Path


VALIDATION_DIR = Path(APPLY_OUTDIR) / "validation"
COHORT_VALIDATION = VALIDATION_DIR / "cohort.json"
RESOURCE_VALIDATION = VALIDATION_DIR / "resources.json"
MODEL_VALIDATION = VALIDATION_DIR / "model.json"
VALIDATION_REPORT = Path(APPLY_OUTDIR) / "validation.json"
VALIDATION_LIBRARY = PYTHON_LIBRARY_DIR / "validation.py"

VALIDATION_VCF = INPUT_VCF_SOURCE
VALIDATION_VCF_INDEX = binding_index(INPUT_VCF_BINDING)
VALIDATION_VCF_INDEX_INPUT = [VALIDATION_VCF_INDEX] if VALIDATION_VCF_INDEX else []
VALIDATION_PED = optional_config_path(INPUT_BINDINGS.get("ped"))
VALIDATION_SAMPLES = optional_config_path(INPUT_BINDINGS.get("samples"))
VALIDATION_FASTA = configured_reference_path("fasta")
VALIDATION_FASTA_INDEX = configured_path(REFERENCE_BINDINGS["fasta"]["fai"])
VALIDATION_FASTA_GZI = optional_config_path(REFERENCE_BINDINGS["fasta"].get("gzi"))
VALIDATION_FASTA_GZI_INPUT = [VALIDATION_FASTA_GZI] if VALIDATION_FASTA_GZI else []

VALIDATION_COVERAGE = optional_reference_path("coverage")
VALIDATION_COVERAGE_INPUT = [VALIDATION_COVERAGE] if VALIDATION_COVERAGE else []
VALIDATION_COVERAGE_INDEX = binding_index(REFERENCE_BINDINGS["coverage"])
VALIDATION_COVERAGE_INDEX_INPUT = [VALIDATION_COVERAGE_INDEX] if VALIDATION_COVERAGE_INDEX else []
VALIDATION_TESTING_UNIT_BED = configured_reference_path("testing_unit_bed")
VALIDATION_TESTING_UNIT_INDEX = binding_index(REFERENCE_BINDINGS["testing_unit_bed"])
VALIDATION_TESTING_UNIT_INDEX_INPUT = [VALIDATION_TESTING_UNIT_INDEX] if VALIDATION_TESTING_UNIT_INDEX else []


def _exclusion_bindings(role):
    return [
        binding
        for metadata, binding in zip(
            config["reference"]["exclusions"], REFERENCE_BINDINGS["exclusions"]
        )
        if metadata["role"] == role
    ]


VALIDATION_VARIANT_EXCLUSION_BINDINGS = _exclusion_bindings("variants_exclude")
VALIDATION_VARIANT_EXCLUSIONS = [binding_data(binding) for binding in VALIDATION_VARIANT_EXCLUSION_BINDINGS]
VALIDATION_VARIANT_EXCLUSION_INDEXES = [binding_index(binding) for binding in VALIDATION_VARIANT_EXCLUSION_BINDINGS]
VALIDATION_PAR_REGION_BINDINGS = _exclusion_bindings("par_regions")
VALIDATION_PAR_REGIONS = [binding_data(binding) for binding in VALIDATION_PAR_REGION_BINDINGS]
VALIDATION_PAR_REGION_INDEXES = [binding_index(binding) for binding in VALIDATION_PAR_REGION_BINDINGS]

VALIDATION_CADD_SNV = binding_data(REFERENCE_BINDINGS["cadd"]["snv"])
VALIDATION_CADD_SNV_INDEX = binding_index(REFERENCE_BINDINGS["cadd"]["snv"])
VALIDATION_CADD_SNV_INDEX_INPUT = [VALIDATION_CADD_SNV_INDEX] if VALIDATION_CADD_SNV_INDEX else []
VALIDATION_CADD_INDEL = binding_data(REFERENCE_BINDINGS["cadd"]["indel"], required=False)
VALIDATION_CADD_INDEL_INPUT = [VALIDATION_CADD_INDEL] if config["model"]["include_indels"] else []
VALIDATION_CADD_INDEL_INDEX = binding_index(REFERENCE_BINDINGS["cadd"]["indel"])
VALIDATION_CADD_INDEL_INDEX_INPUT = (
    [VALIDATION_CADD_INDEL_INDEX]
    if config["model"]["include_indels"] and VALIDATION_CADD_INDEL_INDEX
    else []
)

VALIDATION_LOOKUP_HET = configured_path(REFERENCE_BINDINGS["lookup_tables"]["het"])
VALIDATION_LOOKUP_HOM = configured_path(REFERENCE_BINDINGS["lookup_tables"]["hom"])
VALIDATION_LOOKUP_HEM_CHRX_INPUT = (
    [configured_path(REFERENCE_BINDINGS["lookup_tables"]["hem_chrX"])]
    if config["model"]["hemizygote"]
    else []
)
