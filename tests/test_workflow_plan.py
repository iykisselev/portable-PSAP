"""Unit coverage for immutable workflow mode and binding resolution."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest
import warnings

import yaml


ROOT = Path(__file__).resolve().parents[1]
PYTHON_LIBRARY = ROOT / "workflow" / "lib" / "python"
sys.path.insert(0, str(PYTHON_LIBRARY))

from workflow_plan import (  # noqa: E402
    GENERATED_BUNDLE_CONTRACT,
    WorkflowPlanError,
    build_workflow_plan,
    normalize_run_config,
)


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


class WorkflowPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.science = load_yaml(ROOT / "config" / "config.yaml")

    def manifest(self, relative_path: str) -> dict:
        return merge(self.science, load_yaml(ROOT / relative_path))

    def test_apply_plan_keeps_external_bindings_effective(self) -> None:
        plan = build_workflow_plan(self.manifest("tests/smoke/run.yaml"))

        self.assertEqual("apply", plan.mode)
        self.assertIsNone(plan.generation)
        self.assertIsNotNone(plan.apply)
        assert plan.apply is not None
        self.assertEqual("external", plan.apply.lookup_source.source_type)
        self.assertEqual(
            plan.apply.configured_reference,
            plan.apply.effective_reference,
        )
        self.assertEqual(
            "tests/smoke/data/lookup_het.tsv",
            plan.apply.lookup_source.tables["het"],
        )
        with self.assertRaises(TypeError):
            plan.apply.lookup_source.tables["het"] = "other.tsv"  # type: ignore[index]

    def test_full_plan_derives_generated_lookup_and_effective_scope(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        original = copy.deepcopy(manifest)
        plan = build_workflow_plan(manifest)

        self.assertEqual("full", plan.mode)
        self.assertIsNotNone(plan.generation)
        self.assertIsNotNone(plan.apply)
        assert plan.generation is not None
        assert plan.apply is not None
        self.assertEqual(GENERATED_BUNDLE_CONTRACT, plan.apply.lookup_source.contract)
        self.assertEqual(
            str(plan.generation.bundle_dir / "lookup_hem_chrX.tsv.gz"),
            plan.apply.lookup_source.tables["hem_chrX"],
        )
        self.assertEqual(
            "tests/null_distribution/coverage.bed",
            plan.apply.configured_reference.coverage.data,
        )
        self.assertEqual(
            str(plan.generation.bundle_dir / "calibration_regions.bed.gz"),
            plan.apply.effective_reference.coverage.data,
        )
        self.assertEqual(
            "tests/null_distribution/par_regions.bed",
            plan.generation.par_regions.data,
        )
        self.assertEqual(manifest, original)

    def test_full_plan_disables_hem_lookup_when_the_model_is_disabled(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["model"]["hemizygote"] = False
        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        self.assertNotIn("hem_chrX", plan.apply.lookup_source.tables)

    def test_generate_plan_has_no_apply_inputs_or_external_lookup(self) -> None:
        plan = build_workflow_plan(
            self.manifest("tests/null_distribution/generate.run.yaml")
        )

        self.assertEqual("generate", plan.mode)
        self.assertIsNone(plan.apply)
        self.assertIsNotNone(plan.generation)
        assert plan.generation is not None
        self.assertEqual("AF_CUSTOM", plan.generation.af_source["overall"])
        self.assertEqual(5_000_000, plan.generation.partition_size_bp)

    def test_mutation_rate_policy_resolves_its_declared_binding(self) -> None:
        plan = build_workflow_plan(
            self.manifest("tests/null_distribution/mutation.run.yaml")
        )

        assert plan.generation is not None
        self.assertEqual("mutation_rate", plan.generation.unobserved_policy)
        self.assertIsNotNone(plan.generation.mutation_rate)
        assert plan.generation.mutation_rate is not None
        self.assertEqual(
            "tests/null_distribution/mutation_rates.tsv",
            plan.generation.mutation_rate.data,
        )

    def test_optional_apply_coverage_remains_absent(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["bindings"]["reference"]["coverage"] = {
            "data": None,
            "index": None,
        }
        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        self.assertIsNone(plan.apply.effective_reference.coverage.data)
        self.assertIsNone(plan.apply.lookup_source.calibration_mask)

    def test_full_plan_requires_one_shared_par_role(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["reference"]["exclusions"][1]["role"] = "variants_exclude"

        with self.assertRaisesRegex(WorkflowPlanError, "exactly one.*par_regions"):
            build_workflow_plan(manifest)

    def test_legacy_normalization_warns_without_changing_the_source_mapping(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        del manifest["workflow"]
        original = copy.deepcopy(manifest)

        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            normalized = normalize_run_config(manifest)

        self.assertEqual({"mode": "apply"}, normalized["workflow"])
        self.assertNotIn("workflow", manifest)
        self.assertEqual(original, manifest)
        self.assertEqual(DeprecationWarning, recorded[0].category)


if __name__ == "__main__":
    unittest.main()
