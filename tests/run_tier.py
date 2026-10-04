#!/usr/bin/env python3
"""Run one documented Portable PSAP test tier from the repository root."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]

# Keep the fast loop free of Snakemake subprocesses.  The workflow, E2E, and
# resilience tiers deliberately grow in cost and are selected explicitly.
TIERS: dict[str, tuple[str, ...]] = {
    "unit": (
        "test_annotation_qc.py",
        "test_apply_resources.py",
        "test_bgzf.py",
        "test_bundle_finalization.py",
        "test_combine_exclusions.py",
        "test_execution_policy.py",
        "test_generation_manifests.py",
        "test_generation_partition.py",
        "test_generation_prep.py",
        "test_generation_probability.py",
        "test_generation_resource_policy.py",
        "test_generation_topology.py",
        "test_genomic_contracts.py",
        "test_io_contracts.py",
        "test_null_distribution_oracle.py",
        "test_r_scripts.py",
        "test_sample_metadata_characterization.py",
        "test_sample_manifest.py",
        "test_validate_inputs.py",
        "test_validation_reporting.py",
    ),
    "schema": ("test_workflow_modes.py", "test_workflow_plan.py", "test_stage8_contracts.py"),
    "workflow": ("test_workflow_dags.py", "test_workflow_lint.py", "test_run_examples.py"),
    "e2e": ("test_apply_workflow.py", "test_generate_workflow.py", "test_full_workflow.py"),
    "resilience": ("test_bundle_failure_workflow.py", "test_generation_recovery.py"),
}


def run_pattern(pattern: str) -> int:
    """Run one test file through discovery so existing import paths remain valid."""
    result = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", pattern],
        cwd=ROOT,
        check=False,
    )
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tier", choices=(*TIERS, "all"))
    args = parser.parse_args()
    if args.tier == "all":
        return run_pattern("test_*.py")
    for pattern in TIERS[args.tier]:
        if run_pattern(pattern):
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
