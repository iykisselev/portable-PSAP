"""Workflow-plan coverage for the partitioned-null configuration contract."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from workflow_plan import WorkflowPlanError, build_workflow_plan  # noqa: E402


def merged_generate_config() -> dict:
    base = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    fixture = yaml.safe_load(
        (ROOT / "tests" / "null_distribution" / "generate.run.yaml").read_text()
    )
    for key, value in fixture.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            base[key].update(value)
        else:
            base[key] = value
    return base


class WorkflowPlanTests(unittest.TestCase):
    def test_generation_preserves_configured_partition_order(self) -> None:
        plan = build_workflow_plan(merged_generate_config())
        self.assertEqual(
            ("cds",),
            tuple(
                partition.identifier
                for partition in plan.generation.testing_unit_partitions
            ),
        )

    def test_partition_ids_reject_filename_ambiguous_dots(self) -> None:
        config = merged_generate_config()
        config["generation"]["testing_unit_partitions"][0]["id"] = "cds.part"
        with self.assertRaisesRegex(WorkflowPlanError, "partition"):
            build_workflow_plan(config)

    def test_apply_rejects_removed_external_lookup_form(self) -> None:
        config = merged_generate_config()
        config["workflow"] = {"mode": "apply"}
        config["bindings"] = {"reference": {"lookup_source": {"type": "external"}}}
        config.pop("generation", None)
        with self.assertRaisesRegex(WorkflowPlanError, "bundle"):
            build_workflow_plan(config)

    def test_duplicate_partition_identity_is_rejected(self) -> None:
        config = merged_generate_config()
        duplicate = deepcopy(config["generation"]["testing_unit_partitions"][0])
        config["generation"]["testing_unit_partitions"].append(duplicate)
        with self.assertRaisesRegex(WorkflowPlanError, "duplicate"):
            build_workflow_plan(config)


if __name__ == "__main__":
    unittest.main()
