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
