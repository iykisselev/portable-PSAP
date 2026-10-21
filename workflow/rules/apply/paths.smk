"""Shared paths for apply rules.

``Snakefile`` creates ``WORKFLOW_PLAN`` and ``WORKFLOW_PYTHON_DIR`` before it
loads this file. The constants below give every apply rule the same output
layout.
"""

from pathlib import Path

if WORKFLOW_PLAN.apply is None:
    raise ValueError("Apply rules require an apply or full WorkflowPlan")

APPLY_PLAN = WORKFLOW_PLAN.apply
APPLY_OUTPUT_DIR = APPLY_PLAN.output_dir
APPLY_LOG_DIR = APPLY_PLAN.log_dir
APPLY_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "apply"

ANNOTATION_DIR = APPLY_OUTPUT_DIR / "annotations"
SAMPLE_MANIFEST_PATH = APPLY_OUTPUT_DIR / "sample_manifest.tsv"
TESTING_UNIT = WORKFLOW_PLAN.model.unit
RAW_DIR = APPLY_OUTPUT_DIR / "raw"
ANNOTATED_DIR = APPLY_OUTPUT_DIR / "annotated"
SCORED_DIR = APPLY_OUTPUT_DIR / "scored"
CANDIDATE_DIR = APPLY_OUTPUT_DIR / "candidates"
COHORT_DIR = APPLY_OUTPUT_DIR / "cohort"
PREPARED_DIR = APPLY_OUTPUT_DIR / "prepared"
PREPARED_LOOKUP_DIR = PREPARED_DIR / "lookups"
BENCHMARK_DIR = APPLY_OUTPUT_DIR / "benchmarks"

COHORT_OUTPUTS = {
    name: str(COHORT_DIR / f"{name}.tsv.gz")
    for name in ("cohort_candidates", "cohort_summary", "control_variants")
}
PARTITION_SUMMARIES = {
    partition_id: str(COHORT_DIR / f"cohort_summary.{partition_id}.tsv.gz")
    for partition_id in APPLY_PLAN.lookup_source.partition_ids
}
