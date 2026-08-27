from pathlib import Path
import re
import math
import shlex

import yaml

import pandas as pd


OUTPUT_ROOT = Path(config["outdir"])
APPLY_OUTDIR = OUTPUT_ROOT / "apply" / config["run"]["apply_prefix"]
LOG_DIR = OUTPUT_ROOT / "logs"
REFERENCE_BINDINGS = config["bindings"]["reference"]
INPUT_BINDINGS = config["inputs"]
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts"
R_LIB_DIR = Path(workflow.basedir) / "lib" / "r"
CONDA_ENV_DIR = Path(workflow.basedir) / "envs"
DEPLOYMENT_DIR = Path(workflow.basedir)


def _load_deployment_map(filename):
    with (DEPLOYMENT_DIR / filename).open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


ENVMODULES = _load_deployment_map("envmodules.yaml")
CONTAINERS = _load_deployment_map("containers.yaml")


def mod(*tools):
    """Expand logical tool names into an ordered, de-duplicated module list."""
    modules = []
    for tool in tools:
        for module in ENVMODULES[tool]:
            if module not in modules:
                modules.append(module)
    return modules


def sif(tool):
    """Return the logical tool's pinned Apptainer/Singularity image reference."""
    return CONTAINERS[tool]


def _input_size_mb(input_files):
    size_mb = 0.0
    for item in input_files:
        try:
            size_mb += Path(str(item)).stat().st_size / (1024 * 1024)
        except OSError:
            pass
    return size_mb


def scaled_mem_mb(wildcards, input, attempt, baseline_mb, mb_per_gb=0):
    """Scale expensive jobs by declared input size and retry attempt."""
    return int((baseline_mb + math.ceil(_input_size_mb(input) / 1024) * mb_per_gb) * attempt)


def scaled_runtime_minutes(wildcards, input, attempt, baseline_minutes, minutes_per_gb=0):
    """Scale runtime with the same portable, scheduler-independent estimate."""
    return int((baseline_minutes + math.ceil(_input_size_mb(input) / 1024) * minutes_per_gb) * attempt)


def configured_path(value):
    """Expand a required user-configured path without resolving relative paths."""
    path = str(value).strip()
    if not path or path.lower() in {"na", "none", "null", ".", "-"}:
        raise ValueError(f"Required configured path is missing: {value!r}")
    return str(Path(path).expanduser())


def optional_config_path(value):
    if value is None:
        return ""
    path = str(value).strip()
    if path.lower() in {"", "na", "none", "null", ".", "-"}:
        return ""
    return str(Path(path).expanduser())


def configured_reference_path(name):
    """Return one required physical binding for a logical reference resource."""
    return configured_path(REFERENCE_BINDINGS[name]["data"])


def optional_reference_path(name):
    """Return one optional physical binding, represented as YAML null when absent."""
    binding = REFERENCE_BINDINGS.get(name)
    return optional_config_path(binding["data"] if binding else None)


def binding_data(binding, required=True):
    value = binding.get("data")
    return configured_path(value) if required else optional_config_path(value)


def binding_index(binding):
    """Return only an index explicitly supplied by the run binding."""
    return optional_config_path(binding.get("index"))

R_RUNTIME_FILES = [
    str(APPLY_SCRIPT_DIR / "bootstrap.R"),
    str(R_LIB_DIR / "common.R"),
    str(R_LIB_DIR / "table_contracts.R"),
    str(R_LIB_DIR / "variant_identity.R"),
    str(R_LIB_DIR / "disease_models.R"),
    str(R_LIB_DIR / "genotypes.R"),
    str(R_LIB_DIR / "lookup_contracts.R"),
]

ANNOTATION_DIR = APPLY_OUTDIR / "annotations"
SAMPLE_MANIFEST_PATH = APPLY_OUTDIR / "sample_manifest.tsv"

SAMPLE_MANIFEST_COLUMNS = ["family_id", "pid", "sex", "phenotype"]
SAMPLE_PATTERN = r"[^/]+"
TESTING_UNIT = str(config["model"]["unit"])
TESTING_UNIT_PATTERN = re.escape(TESTING_UNIT)
TESTING_UNIT_TABLE = ANNOTATION_DIR / f"{TESTING_UNIT}.tsv.gz"

RAW_DIR = APPLY_OUTDIR / "raw"
STANDARDIZED_DIR = APPLY_OUTDIR / "standardized"
SEX_CLASSIFIED_DIR = APPLY_OUTDIR / "sex_classified"
ANNOTATED_DIR = APPLY_OUTDIR / "annotated"
SCORED_DIR = APPLY_OUTDIR / "scored"
CANDIDATE_DIR = APPLY_OUTDIR / "candidates"
COHORT_DIR = APPLY_OUTDIR / "cohort"
PREPARED_DIR = APPLY_OUTDIR / "prepared"
BENCHMARK_DIR = APPLY_OUTDIR / "benchmarks"

# Staging is selected solely from validated run metadata.  No filesystem probe may
# add, remove, or redirect a staging job while the workflow is being parsed.
INPUT_VCF_BINDING = INPUT_BINDINGS["vcf"]
INPUT_VCF_SOURCE = binding_data(INPUT_VCF_BINDING)
INPUT_VCF_SOURCE_INDEX = binding_index(INPUT_VCF_BINDING)
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


def table_resource_chromosome_map(wildcards, input):
    if not TABLE_RESOURCE_SPECS[wildcards.resource]["chromosome_map"]:
        return ""
    return f"--chromosome-map-vcf {shlex.quote(str(input.cohort[0]))}"


def table_resource_tsv_layout(wildcards):
    layout = TABLE_RESOURCE_SPECS[wildcards.resource]["tsv_layout"]
    return f"--tsv-layout {layout}" if layout else ""


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


def vcf_resource_source_index(wildcards):
    return [INPUT_VCF_SOURCE_INDEX] if wildcards.resource == "cohort" and INPUT_VCF_SOURCE_INDEX else []


def vcf_resource_validation_marker(wildcards):
    return [validation_marker(wildcards)] if wildcards.resource == "cohort" else []


def validation_marker(wildcards):
    return checkpoints.require_valid_inputs.get().output.ok


def validated_manifest_path(wildcards):
    checkpoints.require_valid_inputs.get()
    return str(SAMPLE_MANIFEST_PATH)


def read_sample_manifest(path):
    manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(SAMPLE_MANIFEST_COLUMNS) - set(manifest.columns)
    if missing:
        raise ValueError(
            f"Sample manifest {path!r} has columns {manifest.columns.tolist()!r}; "
            "missing: " + ", ".join(sorted(missing))
        )

    manifest = manifest.loc[:, SAMPLE_MANIFEST_COLUMNS]
    if manifest["pid"].duplicated().any():
        raise ValueError("Sample manifest contains duplicate pid values")
    return manifest


def validated_sample_manifest(wildcards):
    return read_sample_manifest(validated_manifest_path(wildcards))


def sample_metadata(path, pid, expected_phenotype=None):
    manifest = read_sample_manifest(path)
    rows = manifest.loc[manifest["pid"] == str(pid)]
    if len(rows) != 1:
        raise ValueError(f"Sample {pid!r} is not uniquely present in the manifest")

    metadata = rows.iloc[0]
    if (
        expected_phenotype is not None
        and metadata["phenotype"] != expected_phenotype
    ):
        raise ValueError(
            f"Sample {pid!r} has phenotype {metadata['phenotype']!r}; "
            f"expected {expected_phenotype!r}"
        )
    return metadata


def sample_sex(wildcards):
    return sample_metadata(validated_manifest_path(wildcards), wildcards.sample)["sex"]


def affected_sample_sex(wildcards):
    return sample_metadata(
        validated_manifest_path(wildcards),
        wildcards.sample,
        expected_phenotype="affected",
    )["sex"]


def sample_ids_by_phenotype(wildcards, phenotype):
    manifest = validated_sample_manifest(wildcards)
    return manifest.loc[manifest["phenotype"] == phenotype, "pid"].tolist()


def control_annotated_tables(wildcards):
    return expand(
        str(ANNOTATED_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_by_phenotype(wildcards, "unaffected"),
        testing_unit=[TESTING_UNIT],
    )


def control_psap_tables(wildcards):
    return expand(
        str(SCORED_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_by_phenotype(wildcards, "unaffected"),
        testing_unit=[TESTING_UNIT],
    )


def affected_candidate_tables(wildcards):
    return expand(
        str(
            CANDIDATE_DIR / "{sample}.{testing_unit}.tsv.gz"
        ),
        sample=sample_ids_by_phenotype(wildcards, "affected"),
        testing_unit=[TESTING_UNIT],
    )


def control_pids(wildcards):
    manifest = validated_sample_manifest(wildcards)
    return manifest.loc[manifest["phenotype"] == "unaffected", "pid"].tolist()


def control_sexes(wildcards):
    manifest = validated_sample_manifest(wildcards)
    return manifest.loc[manifest["phenotype"] == "unaffected", "sex"].tolist()
