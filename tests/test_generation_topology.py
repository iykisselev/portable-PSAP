"""Unit contracts for static generation paths with no genomic input reads."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_topology import build_generation_execution_layout  # noqa: E402
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
    def test_execution_layout_derives_bins_and_roles_from_the_model(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(ROOT / "tests" / "null_distribution" / "generate.run.yaml"),
        )
        plan = build_workflow_plan(configuration)
        assert plan.generation is not None

        layout = build_generation_execution_layout(
            plan.generation,
            plan.model,
            plan.generation_settings,
        )
        self.assertEqual(1400, layout.bin_count)
        self.assertEqual(("het", "hom", "hem_chrX"), layout.roles)
        self.assertEqual(
            str(plan.generation.bundle_dir / "lookup.cds.DOM-het.tsv.gz"),
            layout.partition_lookup_paths[("cds", "DOM-het")],
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
            str(
                plan.generation.scalable_dir
                / "windows"
                / "1_000000000001_000000005000.allele_states.parquet"
            ),
            context.window_artifact(
                "1_000000000001_000000005000",
                "allele_states",
            ),
        )
        self.assertEqual(
            str(plan.generation.bundle_dir / "lookup_het.tsv.gz"),
            context.bundle_outputs()["lookup_het"],
        )
        self.assertEqual(
            str(plan.generation.bundle_dir / "lookup_hem_chrX.tsv.gz"),
            context.bundle_outputs()["lookup_hem_chrX"],
        )
        self.assertNotIn("validation", context.bundle_outputs())
        self.assertEqual(
            f"population_vcf:{plan.generation.af_source['overall']}",
            context.bundle_resource_identities["population_vcf"]["identity"],
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

    def test_rule_context_keeps_optional_hem_output_named_when_disabled(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(
                ROOT / "tests" / "null_distribution" / "generate_no_hem.run.yaml"
            ),
        )
        context = GenerationRuleContext.from_workflow_plan(
            build_workflow_plan(configuration)
        )

        self.assertFalse(context.model.hemizygote)
        self.assertEqual([], context.bundle_outputs()["lookup_hem_chrX"])
        self.assertNotIn("marker", context.bundle_outputs())


if __name__ == "__main__":
    unittest.main()
