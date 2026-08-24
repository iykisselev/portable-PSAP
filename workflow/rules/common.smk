from pathlib import Path
import math
import re

import pandas as pd
import yaml


OUTPUT_ROOT = Path(config["outdir"])
APPLY_OUTDIR = OUTPUT_ROOT / "apply" / config["run"]["apply_prefix"]
LOG_DIR = OUTPUT_ROOT / "logs"
REFERENCE_BINDINGS = config["bindings"]["reference"]
INPUT_BINDINGS = config["inputs"]
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts"
PYTHON_LIBRARY_DIR = Path(workflow.basedir) / "lib" / "python"
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
    """Return an optional physical binding, represented as an empty value when absent."""
    binding = REFERENCE_BINDINGS.get(name)
    return optional_config_path(binding["data"] if binding else None)


def binding_data(binding, required=True):
    value = binding.get("data")
    return configured_path(value) if required else optional_config_path(value)


def binding_index(binding):
    """Return only an index explicitly supplied by the run binding."""
    return optional_config_path(binding.get("index"))


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
        str(CANDIDATE_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_by_phenotype(wildcards, "affected"),
        testing_unit=[TESTING_UNIT],
    )


def control_pids(wildcards):
    manifest = validated_sample_manifest(wildcards)
    return manifest.loc[manifest["phenotype"] == "unaffected", "pid"].tolist()


def control_sexes(wildcards):
    manifest = validated_sample_manifest(wildcards)
    return manifest.loc[manifest["phenotype"] == "unaffected", "sex"].tolist()
