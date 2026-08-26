from pathlib import Path
import re

import pandas as pd


BGZF_EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)

OUTPUT_ROOT = Path(config["outdir"])
APPLY_OUTDIR = OUTPUT_ROOT / "apply" / config["run"]["apply_prefix"]
LOG_DIR = OUTPUT_ROOT / "logs"
REFERENCE_BINDINGS = config["bindings"]["reference"]
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "apply"
R_LIB_DIR = Path(workflow.basedir) / "lib" / "r"
CONDA_ENV_DIR = Path(workflow.basedir) / "envs"
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


def tabix_index(path):
    path = str(path)
    tbi = Path(f"{path}.tbi")
    return str(tbi if tbi.is_file() else Path(f"{path}.csi"))


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
    return configured_path(REFERENCE_BINDINGS[name])


def optional_reference_path(name):
    """Return one optional physical binding, represented as YAML null when absent."""
    return optional_config_path(REFERENCE_BINDINGS.get(name))


def comma_separated_config_paths(value):
    """Parse a nullable, comma-separated config value into ordered paths."""
    normalized = optional_config_path(value)
    if not normalized:
        return []
    return [
        path
        for item in normalized.split(",")
        if (path := optional_config_path(item))
    ]


def existing_genomic_indexes(paths):
    """Return existing tabix/CSI sidecars inspected by validation rules."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    return [
        str(index)
        for path in paths
        for suffix in (".tbi", ".csi")
        if (index := Path(f"{path}{suffix}")).is_file()
    ]


def genomic_resource_kind(path):
    """Classify supported genomic resources independently of compression."""
    name = str(path).lower()
    for suffix in (".gz", ".bgz", ".bgzf"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith(".vcf"):
        return "vcf"
    if name.endswith(".bed"):
        return "bed"
    if name.endswith((".tsv", ".txt")):
        return "tsv"
    return "unknown"


def is_bgzip_file(path):
    """Return whether path starts with a valid BGZF gzip extra block."""
    try:
        with Path(path).open("rb") as handle:
            header = handle.read(18)
            if len(header) < 18 or header[:2] != b"\x1f\x8b" or not header[3] & 4:
                return False
            extra_length = int.from_bytes(header[10:12], "little")
            extra = header[12:] + handle.read(max(0, extra_length - 6))
    except OSError:
        return False

    offset = 0
    while offset + 4 <= len(extra):
        subfield_length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        if extra[offset : offset + 2] == b"BC" and subfield_length == 2:
            return True
        offset += 4 + subfield_length
    return False


def has_genomic_index(path):
    return any(Path(f"{path}{suffix}").is_file() for suffix in (".tbi", ".csi"))


def has_bgzf_eof(path):
    try:
        with Path(path).open("rb") as handle:
            handle.seek(-len(BGZF_EOF_BLOCK), 2)
            return handle.read() == BGZF_EOF_BLOCK
    except OSError:
        return False


def genomic_resource_is_ready(path):
    """Resources already BGZF-compressed and indexed can be used in place."""
    return (
        Path(path).is_file()
        and is_bgzip_file(path)
        and has_bgzf_eof(path)
        and has_genomic_index(path)
    )


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
