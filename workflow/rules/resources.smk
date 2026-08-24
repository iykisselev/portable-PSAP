# Build the resource registry from scientific identities and one run binding.
# These helpers deliberately select jobs only from configuration, never probes.
INPUT_VCF_BINDING = INPUT_BINDINGS["vcf"]
INPUT_VCF_SOURCE = binding_data(INPUT_VCF_BINDING)
PREPARED_INPUT_VCF = str(PREPARED_DIR / "cohort.vcf.gz")
WORKFLOW_INPUT_VCF = PREPARED_INPUT_VCF
WORKFLOW_INPUT_VCF_INDEX = f"{PREPARED_INPUT_VCF}.tbi"

PREPARED_TABLE_DIR = PREPARED_DIR / "tables"
TABLE_RESOURCE_SPECS = {}


def add_table_resource(resource, binding, kind, *, chromosome_map=False, tsv_layout=None):
    data = binding_data(binding, required=False)
    if not data:
        return
    TABLE_RESOURCE_SPECS[resource] = {
        "source": data,
        "source_index": binding_index(binding),
        "kind": kind,
        "chromosome_map": chromosome_map,
        "tsv_layout": tsv_layout,
    }


add_table_resource(
    "coverage",
    REFERENCE_BINDINGS["coverage"],
    config["reference"]["coverage"]["kind"],
    chromosome_map=True,
)
add_table_resource(
    "testing_units", REFERENCE_BINDINGS["testing_unit_bed"], "bed", chromosome_map=True
)
add_table_resource(
    "cadd_snv", REFERENCE_BINDINGS["cadd"]["snv"], "tsv", chromosome_map=True, tsv_layout="position"
)
if config["model"]["include_indels"]:
    add_table_resource(
        "cadd_indel", REFERENCE_BINDINGS["cadd"]["indel"], "tsv", chromosome_map=True, tsv_layout="position"
    )

EXCLUSION_IDS = []
EXCLUSION_SOURCE_BY_ID = {}
EXCLUSION_KIND_BY_ID = {}
if len(config["reference"]["exclusions"]) != len(REFERENCE_BINDINGS["exclusions"]):
    raise ValueError(
        "reference.exclusions and bindings.reference.exclusions must have the same length"
    )
for position, (metadata, binding) in enumerate(
    zip(config["reference"]["exclusions"], REFERENCE_BINDINGS["exclusions"])
):
    item = f"exclusion_{position:03d}"
    add_table_resource(item, binding, metadata["kind"])
    if item in TABLE_RESOURCE_SPECS:
        EXCLUSION_IDS.append(item)
        EXCLUSION_SOURCE_BY_ID[item] = TABLE_RESOURCE_SPECS[item]["source"]
        EXCLUSION_KIND_BY_ID[item] = metadata["kind"]

TABLE_RESOURCE_IDS = list(TABLE_RESOURCE_SPECS)
TABLE_RESOURCE_PATTERN = f"(?:{'|'.join(map(re.escape, TABLE_RESOURCE_IDS))})"
TABLE_RESOURCE_DATA = {
    resource: str(PREPARED_TABLE_DIR / f"{resource}.gz")
    for resource in TABLE_RESOURCE_IDS
}
TABLE_RESOURCE_INDEX = {
    resource: f"{path}.tbi" for resource, path in TABLE_RESOURCE_DATA.items()
}
TABLE_RESOURCE_DEFAULTS = {
    "coverage": (2048, 60),
    "testing_units": (2048, 60),
    "cadd_snv": (2048, 360),
    "cadd_indel": (2048, 360),
}

WORKFLOW_COVERAGE = TABLE_RESOURCE_DATA.get("coverage", "")
WORKFLOW_COVERAGE_INDEX = TABLE_RESOURCE_INDEX.get("coverage", "")
WORKFLOW_TESTING_UNIT = TABLE_RESOURCE_DATA["testing_units"]
WORKFLOW_TESTING_UNIT_INDEX = TABLE_RESOURCE_INDEX["testing_units"]
WORKFLOW_CADD_SNV = TABLE_RESOURCE_DATA["cadd_snv"]
WORKFLOW_CADD_SNV_INDEX = TABLE_RESOURCE_INDEX["cadd_snv"]
WORKFLOW_CADD_INDEL = TABLE_RESOURCE_DATA.get("cadd_indel", "")
WORKFLOW_CADD_INDEL_INDEX = TABLE_RESOURCE_INDEX.get("cadd_indel", "")
WORKFLOW_EXCLUSIONS = [TABLE_RESOURCE_DATA[item] for item in EXCLUSION_IDS]
WORKFLOW_EXCLUSION_INDEXES = [TABLE_RESOURCE_INDEX[item] for item in EXCLUSION_IDS]
EXCLUSION_SOURCES = [EXCLUSION_SOURCE_BY_ID[item] for item in EXCLUSION_IDS]

COMBINED_EXCLUSIONS_VCF = str(PREPARED_DIR / "combined_exclusions.vcf")
COMBINED_EXCLUSIONS = str(PREPARED_DIR / "combined_exclusions.vcf.gz")
COMBINED_EXCLUSIONS_INDEX = f"{COMBINED_EXCLUSIONS}.tbi"
VCF_RESOURCE_IDS = ["cohort", "combined_exclusions"]
VCF_RESOURCE_PATTERN = f"(?:{'|'.join(VCF_RESOURCE_IDS)})"


def table_resource_source(wildcards):
    return TABLE_RESOURCE_SPECS[wildcards.resource]["source"]


def table_resource_source_index(wildcards):
    index = TABLE_RESOURCE_SPECS[wildcards.resource]["source_index"]
    return [index] if index else []


def table_resource_cohort(wildcards):
    return [WORKFLOW_INPUT_VCF] if TABLE_RESOURCE_SPECS[wildcards.resource]["chromosome_map"] else []


def table_resource_cohort_index(wildcards):
    return [WORKFLOW_INPUT_VCF_INDEX] if TABLE_RESOURCE_SPECS[wildcards.resource]["chromosome_map"] else []


def table_resource_mem_mb(wildcards, input, attempt):
    baseline_mb, _ = TABLE_RESOURCE_DEFAULTS.get(wildcards.resource, (4096, 180))
    return scaled_mem_mb(
        wildcards,
        input,
        attempt,
        baseline_mb,
        256 if wildcards.resource.startswith("cadd_") else 0,
    )


def table_resource_runtime_minutes(wildcards, input, attempt):
    _, baseline_minutes = TABLE_RESOURCE_DEFAULTS.get(wildcards.resource, (4096, 180))
    return scaled_runtime_minutes(
        wildcards,
        input,
        attempt,
        baseline_minutes,
        20 if wildcards.resource.startswith("cadd_") else 0,
    )


def vcf_resource_source(wildcards):
    if wildcards.resource == "cohort":
        return INPUT_VCF_SOURCE
    return rules.combine_exclusions.output.vcf


def vcf_resource_validation_marker(wildcards):
    return [validation_marker(wildcards)] if wildcards.resource == "cohort" else []
