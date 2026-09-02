"""Plan-derived apply paths and resource-sizing adapters.

Requires: ``WORKFLOW_PLAN`` and ``PYTHON_LIBRARY_DIR`` from ``Snakefile``.
Exports only apply paths plus the Snakemake-compatible resource callbacks used
by apply rule modules.
"""

from pathlib import Path
import re

from workflow_runtime import scaled_mem_mb as _scaled_mem_mb
from workflow_runtime import scaled_runtime_minutes as _scaled_runtime_minutes
from execution_policy import resource_mem_mb as _resource_mem_mb
from execution_policy import resource_runtime_minutes as _resource_runtime_minutes


if WORKFLOW_PLAN.apply is None:
    raise ValueError("Apply rules require an apply or full WorkflowPlan")

APPLY_PLAN = WORKFLOW_PLAN.apply
APPLY_OUTDIR = APPLY_PLAN.output_dir
APPLY_LOG_DIR = APPLY_PLAN.log_dir
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts"

ANNOTATION_DIR = APPLY_OUTDIR / "annotations"
SAMPLE_MANIFEST_PATH = APPLY_OUTDIR / "sample_manifest.tsv"
TESTING_UNIT = WORKFLOW_PLAN.model.unit
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


def scaled_mem_mb(wildcards, input, attempt, baseline_mb, mb_per_gb=0):
    """Adapt pure resource sizing to Snakemake's callback signature."""
    return _scaled_mem_mb(input, attempt, baseline_mb, mb_per_gb)


def scaled_runtime_minutes(wildcards, input, attempt, baseline_minutes, minutes_per_gb=0):
    """Adapt pure runtime sizing to Snakemake's callback signature."""
    return _scaled_runtime_minutes(input, attempt, baseline_minutes, minutes_per_gb)


def execution_mem_mb(wildcards, input, attempt, policy_name):
    """Adapt a named execution policy to Snakemake's callback signature."""
    return _resource_mem_mb(policy_name, input, attempt)


def execution_runtime_minutes(wildcards, input, attempt, policy_name):
    """Adapt a named execution policy to Snakemake's callback signature."""
    return _resource_runtime_minutes(policy_name, input, attempt)
