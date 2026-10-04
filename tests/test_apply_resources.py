"""Contracts for immutable apply resource-registry construction."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from apply_resources import build_apply_resource_registry  # noqa: E402
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


class ApplyResourceRegistryTests(unittest.TestCase):
    def test_primary_exclusion_is_subtracted_not_used_as_calibration_regions(self):
        configuration = merge(load_yaml(ROOT / "config/config.yaml"), load_yaml(ROOT / "tests/smoke/run.yaml"))
        configuration["bindings"]["reference"]["exclusion"] = {"data": "mask.bed", "index": None}
        registry = build_apply_resource_registry(build_workflow_plan(configuration), Path("output/prepared"))
        self.assertEqual(registry.resource("exclusion").scientific_role, "variants_exclude")
        self.assertIn("exclusion", registry.exclusion_ids)
        self.assertNotIn("calibration_regions", registry.resources)

    def test_registry_uses_one_effective_binding_set_and_is_immutable(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(ROOT / "tests" / "smoke" / "run.yaml"),
        )
        plan = build_workflow_plan(configuration)
        registry = build_apply_resource_registry(plan, Path("output") / "prepared")

        self.assertEqual(
            ("testing_units", "cadd_snv", "cadd_indel", "exclusion_000", "exclusion_001"),
            registry.table_resource_ids,
        )
        self.assertEqual(
            "(?:testing_units|cadd_snv|cadd_indel|exclusion_000|exclusion_001)",
            registry.table_resource_pattern,
        )
        self.assertEqual("output/prepared/cohort.vcf.gz", registry.prepared_input_vcf)
        self.assertEqual(
            "tests/smoke/data/exclude_full_span.bed",
            registry.exclusion_source_by_id["exclusion_000"],
        )
        with self.assertRaises(TypeError):
            registry.table_resources["other"] = registry.table_resources["cadd_snv"]  # type: ignore[index]

        cadd = registry.resource("cadd_snv")
        self.assertEqual("cadd_snv", cadd.scientific_role)
        self.assertTrue(cadd.requires_chromosome_normalization)
        self.assertEqual(cadd.source, "tests/smoke/data/cadd_snv.tsv")
        self.assertEqual(
            cadd.staged_data,
            "output/prepared/tables/cadd_snv.tsv.gz",
        )

    def test_full_registry_distinguishes_bundle_scope_from_shared_inputs(self) -> None:
        configuration = merge(
            load_yaml(ROOT / "config" / "config.yaml"),
            load_yaml(ROOT / "tests" / "null_distribution" / "full.run.yaml"),
        )
        plan = build_workflow_plan(configuration)
        registry = build_apply_resource_registry(plan, Path("output") / "prepared")

        assert plan.apply is not None
        assert plan.generation is not None
        self.assertFalse(registry.resource("fasta").stage)
        self.assertEqual("configured", registry.resource("fasta").origin)
        self.assertEqual("configured", registry.resource("cadd_snv").origin)
        self.assertEqual(
            plan.apply.configured_reference.cadd_snv.data,
            registry.resource("cadd_snv").source,
        )
        self.assertEqual("generated_bundle", registry.resource("testing_units").origin)
        self.assertEqual("generated_bundle", registry.resource("calibration_regions").origin)
        self.assertEqual("generated_bundle", registry.resource("exclusion_000").origin)
        self.assertEqual(
            str(plan.generation.bundle_dir / "eligible_testing_units.bed.gz"),
            registry.resource("testing_units").source,
        )


if __name__ == "__main__":
    unittest.main()
