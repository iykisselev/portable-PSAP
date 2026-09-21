"""Shared paths and resource callbacks for apply rules.

``Snakefile`` creates ``WORKFLOW_PLAN`` and ``WORKFLOW_PYTHON_DIR`` before it
loads this file. The constants below give every apply rule the same output
layout and resource-sizing functions.
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
APPLY_OUTPUT_DIR = APPLY_PLAN.output_dir
APPLY_LOG_DIR = APPLY_PLAN.log_dir
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts"

ANNOTATION_DIR = APPLY_OUTPUT_DIR / "annotations"
SAMPLE_MANIFEST_PATH = APPLY_OUTPUT_DIR / "sample_manifest.tsv"
TESTING_UNIT = WORKFLOW_PLAN.model.unit
TESTING_UNIT_PATTERN = re.escape(TESTING_UNIT)
TESTING_UNIT_TABLE = ANNOTATION_DIR / f"{TESTING_UNIT}.tsv.gz"

RAW_DIR = APPLY_OUTPUT_DIR / "raw"
STANDARDIZED_DIR = APPLY_OUTPUT_DIR / "standardized"
SEX_CLASSIFIED_DIR = APPLY_OUTPUT_DIR / "sex_classified"
ANNOTATED_DIR = APPLY_OUTPUT_DIR / "annotated"
SCORED_DIR = APPLY_OUTPUT_DIR / "scored"
CANDIDATE_DIR = APPLY_OUTPUT_DIR / "candidates"
COHORT_DIR = APPLY_OUTPUT_DIR / "cohort"
PREPARED_DIR = APPLY_OUTPUT_DIR / "prepared"
BENCHMARK_DIR = APPLY_OUTPUT_DIR / "benchmarks"


def scaled_mem_mb(wildcards, input, attempt, baseline_mb, mb_per_gb=0):
    """Use the shared memory calculation with Snakemake's callback arguments."""
    return _scaled_mem_mb(input, attempt, baseline_mb, mb_per_gb)


def scaled_runtime_minutes(wildcards, input, attempt, baseline_minutes, minutes_per_gb=0):
    """Use the shared runtime calculation with Snakemake's callback arguments."""
    return _scaled_runtime_minutes(input, attempt, baseline_minutes, minutes_per_gb)


def execution_mem_mb(wildcards, input, attempt, policy_name):
    """Look up a named memory policy for this Snakemake job."""
    return _resource_mem_mb(policy_name, input, attempt)


def execution_runtime_minutes(wildcards, input, attempt, policy_name):
    """Look up a named runtime policy for this Snakemake job."""
    return _resource_runtime_minutes(policy_name, input, attempt)
