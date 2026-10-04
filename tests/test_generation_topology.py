"""Unit contracts for static generation paths with no genomic input reads."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_topology import (  # noqa: E402
    build_generation_execution_layout,
    build_generation_scope_artifacts,
)
from generation_context import GenerationRuleContext  # noqa: E402
from workflow_plan import build_workflow_plan  # noqa: E402


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def merge(base: dict, override: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


class GenerationTopologyTests(unittest.TestCase):
    def test_static_context_defers_topology_to_scope_outputs(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(ROOT / "tests" / "null_distribution" / "generate.run.yaml"),
        )
        plan = build_workflow_plan(configuration)
        assert plan.generation is not None

        scope = build_generation_scope_artifacts(plan.generation)
        self.assertEqual(
            str(plan.generation.scalable_dir / "prepared" / "partitions.json"),
            scope.artifacts["partitions"],
        )
        self.assertEqual(
            str(plan.generation.scalable_dir / "prepared" / "unit_buckets.json"),
            scope.artifacts["buckets"],
        )

        layout = build_generation_execution_layout(
            plan.generation,
            plan.model,
            plan.generation_settings,
        )
        self.assertEqual(1400, layout.bin_count)
        self.assertEqual(("het", "hom", "hem_chrX"), layout.roles)
        self.assertEqual(
            str(plan.generation.bundle_dir / "lookup_het.tsv.gz"),
            layout.bundle_paths["lookup_het"],
        )

    def test_rule_context_keeps_generation_layout_and_inputs_together(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(ROOT / "tests" / "null_distribution" / "generate.run.yaml"),
        )
        plan = build_workflow_plan(configuration)
        context = GenerationRuleContext.from_workflow_plan(plan)

        self.assertEqual(plan.generation, context.plan)
        self.assertEqual(
            str(plan.generation.scalable_dir / "partitions"),
            str(context.partition_dir),
        )
        self.assertEqual(
            str(plan.generation.scalable_dir / "calculated" / "partition_statistics"),
            str(context.statistics_dir),
        )
        self.assertEqual(
            str(
                plan.generation.scalable_dir
                / "partitions"
                / "1_000000000001_000000005000.allele_states.parquet"
            ),
            context.partition_artifact(
                "1_000000000001_000000005000",
                "allele_states",
            ),
        )
        self.assertEqual(
            str(
                plan.generation.scalable_dir
                / "calculated"
                / "lookup_buckets"
                / "000000.het.tsv"
            ),
            context.bucket_artifact("000000", "lookup_het"),
        )


if __name__ == "__main__":
    unittest.main()
