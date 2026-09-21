"""Inputs shared by the apply-validation rules.

The earlier apply files have already resolved every path. Keeping those values
here lets each validation rule use the same sources as the downstream rules.
"""

from pathlib import Path


VALIDATION_DIR = Path(APPLY_OUTPUT_DIR) / "validation"
COHORT_VALIDATION_REPORT = VALIDATION_DIR / "cohort.json"
RESOURCE_VALIDATION_REPORT = VALIDATION_DIR / "resources.json"
MODEL_VALIDATION_REPORT = VALIDATION_DIR / "model.json"
COMBINED_VALIDATION_REPORT = Path(APPLY_OUTPUT_DIR) / "validation.json"
VALIDATION_LIBRARIES = [
    str(WORKFLOW_PYTHON_DIR / filename)
    for filename in (
        "validation.py",
        "validation_cohort.py",
        "validation_resources.py",
        "validation_model.py",
        "validation_reporting.py",
        "bundle_certificate.py",
        "sample_manifest.py",
    )
]
VALIDATION_REPORT_SCHEMA = str(
    Path(workflow.basedir) / "schemas" / "validation_report.schema.yaml"
)

def required_resource(name):
    return RESOURCE_SPECS[name]


def optional_resource(name):
    return RESOURCE_SPECS.get(name)


VALIDATION_VCF = required_resource("cohort").source
VALIDATION_VCF_INDEX = required_resource("cohort").source_index
VALIDATION_VCF_INDEX_INPUT = [VALIDATION_VCF_INDEX] if VALIDATION_VCF_INDEX else []
VALIDATION_PED = INPUT_BINDINGS.ped or ""
VALIDATION_SAMPLES = INPUT_BINDINGS.samples or ""
VALIDATION_FASTA = required_resource("fasta").source
VALIDATION_FASTA_INDEX = required_resource("fasta").source_index
VALIDATION_FASTA_GZI = (
    required_resource("fasta").companion_indexes[0]
    if required_resource("fasta").companion_indexes
    else ""
)
VALIDATION_FASTA_GZI_INPUT = [VALIDATION_FASTA_GZI] if VALIDATION_FASTA_GZI else []

VALIDATION_COVERAGE_SPEC = optional_resource("coverage")
VALIDATION_COVERAGE = (
    VALIDATION_COVERAGE_SPEC.source if VALIDATION_COVERAGE_SPEC else ""
)
VALIDATION_COVERAGE_INPUT = [VALIDATION_COVERAGE] if VALIDATION_COVERAGE else []
VALIDATION_COVERAGE_INDEX = (
    VALIDATION_COVERAGE_SPEC.source_index if VALIDATION_COVERAGE_SPEC else ""
)
VALIDATION_COVERAGE_INDEX_INPUT = [VALIDATION_COVERAGE_INDEX] if VALIDATION_COVERAGE_INDEX else []
VALIDATION_TESTING_UNIT_BED = required_resource("testing_units").source
VALIDATION_TESTING_UNIT_INDEX = required_resource("testing_units").source_index
VALIDATION_TESTING_UNIT_INDEX_INPUT = [VALIDATION_TESTING_UNIT_INDEX] if VALIDATION_TESTING_UNIT_INDEX else []


def _exclusion_resources(role):
    return [
        spec
        for spec in RESOURCE_SPECS.values()
        if spec.scientific_role == role
    ]


VALIDATION_VARIANT_EXCLUSION_SPECS = _exclusion_resources("variants_exclude")
VALIDATION_VARIANT_EXCLUSIONS = [spec.source for spec in VALIDATION_VARIANT_EXCLUSION_SPECS]
VALIDATION_VARIANT_EXCLUSION_INDEXES = [spec.source_index for spec in VALIDATION_VARIANT_EXCLUSION_SPECS]
VALIDATION_PAR_REGION_SPECS = _exclusion_resources("par_regions")
VALIDATION_PAR_REGIONS = [spec.source for spec in VALIDATION_PAR_REGION_SPECS]
VALIDATION_PAR_REGION_INDEXES = [spec.source_index for spec in VALIDATION_PAR_REGION_SPECS]

VALIDATION_CADD_SNV = required_resource("cadd_snv").source
VALIDATION_CADD_SNV_INDEX = required_resource("cadd_snv").source_index
VALIDATION_CADD_SNV_INDEX_INPUT = [VALIDATION_CADD_SNV_INDEX] if VALIDATION_CADD_SNV_INDEX else []
VALIDATION_CADD_INDEL_SPEC = optional_resource("cadd_indel")
VALIDATION_CADD_INDEL = (
    VALIDATION_CADD_INDEL_SPEC.source if VALIDATION_CADD_INDEL_SPEC else ""
)
VALIDATION_CADD_INDEL_INPUT = [VALIDATION_CADD_INDEL] if WORKFLOW_PLAN.model.include_indels else []
VALIDATION_CADD_INDEL_INDEX = (
    VALIDATION_CADD_INDEL_SPEC.source_index if VALIDATION_CADD_INDEL_SPEC else ""
)
VALIDATION_CADD_INDEL_INDEX_INPUT = (
    [VALIDATION_CADD_INDEL_INDEX]
    if WORKFLOW_PLAN.model.include_indels and VALIDATION_CADD_INDEL_INDEX
    else []
)

VALIDATION_LOOKUP_HET = LOOKUP_SOURCE.tables["het"]
VALIDATION_LOOKUP_HOM = LOOKUP_SOURCE.tables["hom"]
VALIDATION_LOOKUP_HEM_CHRX_INPUT = (
    [LOOKUP_SOURCE.tables["hem_chrX"]]
    if WORKFLOW_PLAN.model.hemizygote
    else []
)
VALIDATION_LOOKUP_METADATA_INPUT = [LOOKUP_SOURCE.metadata] if LOOKUP_SOURCE.metadata else []
VALIDATION_LOOKUP_REPORT_INPUT = [LOOKUP_SOURCE.validation] if LOOKUP_SOURCE.validation else []
VALIDATION_LOOKUP_MARKER_INPUT = [LOOKUP_SOURCE.marker] if LOOKUP_SOURCE.marker else []
VALIDATION_LOOKUP_ELIGIBLE_UNITS_INPUT = [LOOKUP_SOURCE.eligible_units] if LOOKUP_SOURCE.is_bundle else []
VALIDATION_LOOKUP_CALIBRATION_MASK_INPUT = [LOOKUP_SOURCE.calibration_mask] if LOOKUP_SOURCE.is_bundle else []
VALIDATION_LOOKUP_EXCLUDED_ALLELES_INPUT = list(LOOKUP_SOURCE.excluded_alleles) if LOOKUP_SOURCE.is_bundle else []
VALIDATION_LOOKUP_EXCLUDED_UNITS_INPUT = [LOOKUP_SOURCE.excluded_units] if LOOKUP_SOURCE.is_bundle else []
