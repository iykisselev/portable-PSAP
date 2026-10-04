"""Configuration contracts for staged null-distribution workflow modes."""

from __future__ import annotations

import copy
import unittest
from pathlib import Path

import jsonschema
import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_SCHEMA = ROOT / "workflow" / "schemas" / "run.schema.yaml"
SCIENCE_SCHEMA = ROOT / "workflow" / "schemas" / "science.schema.yaml"


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle)


def merge(base: dict, override: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def validate_run_config(config: dict) -> list[jsonschema.ValidationError]:
    schema = load_yaml(RUN_SCHEMA)
    resolver = jsonschema.RefResolver(base_uri=RUN_SCHEMA.resolve().as_uri(), referrer=schema)
    validator = jsonschema.Draft202012Validator(schema, resolver=resolver)
    return sorted(validator.iter_errors(config), key=lambda error: list(error.path))


def validate_science_config(config: dict) -> list[jsonschema.ValidationError]:
    schema = load_yaml(SCIENCE_SCHEMA)
    validator = jsonschema.Draft202012Validator(schema)
    return sorted(validator.iter_errors(config), key=lambda error: list(error.path))


class WorkflowModeSchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.science = load_yaml(ROOT / "config" / "config.yaml")

    def merged_manifest(self, relative_path: str) -> dict:
        return merge(self.science, load_yaml(ROOT / relative_path))

    def test_legacy_apply_manifest_is_valid_after_normalization(self) -> None:
        config = self.merged_manifest("tests/smoke/run.yaml")
        config["workflow"] = {"mode": "apply"}
        self.assertEqual([], validate_run_config(config))

    def test_generate_manifest_needs_no_cohort_or_external_lookup(self) -> None:
        config = self.merged_manifest("tests/null_distribution/generate.run.yaml")
        self.assertNotIn("inputs", config)
        self.assertNotIn("bindings", config)
        self.assertEqual([], validate_run_config(config))

    def test_science_schema_rejects_removed_reference_metadata(self) -> None:
        config = copy.deepcopy(self.science)
        config["reference"]["lookup_tables"] = {}
        self.assertTrue(validate_science_config(config))

        config = copy.deepcopy(self.science)
        config["reference"]["exclusion"]["optional"] = True
        self.assertTrue(validate_science_config(config))

    def test_full_manifest_uses_generation_inputs_without_external_lookup(self) -> None:
        config = self.merged_manifest("tests/null_distribution/full.run.yaml")
        self.assertNotIn("lookup_source", config["bindings"]["reference"])
        self.assertEqual(
            {"population_vcf", "af_source", "unobserved_policy", "execution"},
            set(config["generation"]),
        )
        self.assertEqual([], validate_run_config(config))

    def test_full_rejects_a_second_shared_resource_binding_set(self) -> None:
        config = self.merged_manifest("tests/null_distribution/full.run.yaml")
        config["generation"]["score_snv"] = config["bindings"]["reference"]["cadd"]["snv"]
        self.assertTrue(validate_run_config(config))

    def test_full_allows_no_exclusion_in_the_shared_reference_bindings(self) -> None:
        config = self.merged_manifest("tests/null_distribution/full.run.yaml")
        config["bindings"]["reference"]["exclusion"] = {"data": None, "index": None}
        self.assertEqual([], validate_run_config(config))

    def test_sex_specific_hem_generation_requires_male_af_tag(self) -> None:
        config = self.merged_manifest("tests/null_distribution/generate.run.yaml")
        del config["generation"]["af_source"]["male"]
        self.assertTrue(validate_run_config(config))

    def test_generate_rejects_apply_resources(self) -> None:
        config = self.merged_manifest("tests/null_distribution/generate.run.yaml")
        config["inputs"] = {"vcf": {"data": "cohort.vcf", "index": None}, "ped": "cohort.ped", "samples": None}
        errors = validate_run_config(config)
        self.assertTrue(errors)
        self.assertIn(
            "not valid under any of the given schemas",
            "\n".join(error.message for error in errors),
        )

    def test_full_rejects_external_lookup_source(self) -> None:
        config = self.merged_manifest("tests/null_distribution/full.run.yaml")
        config["bindings"]["reference"]["lookup_source"] = {
            "type": "external",
            "contract": "external_legacy_v0",
            "tables": {"het": "het.tsv", "hom": "hom.tsv"},
        }
        self.assertTrue(validate_run_config(config))

    def test_apply_without_hemizygote_does_not_require_hem_lookup(self) -> None:
        config = self.merged_manifest("tests/smoke/run.yaml")
        config["workflow"] = {"mode": "apply"}
        config["model"]["hemizygote"] = False
        del config["bindings"]["reference"]["lookup_source"]["tables"]["hem_chrX"]
        self.assertEqual([], validate_run_config(config))

    def test_bundle_lookup_source_requires_every_certified_scope_component(self) -> None:
        config = self.merged_manifest("tests/smoke/run.yaml")
        config["workflow"] = {"mode": "apply"}
        config["model"]["hemizygote"] = False
        config["model"]["include_indels"] = False
        source = {
            "type": "bundle",
            "contract": "generated_bundle_v1",
            "tables": {"het": "lookup_het.tsv.gz", "hom": "lookup_hom.tsv.gz"},
            "metadata": "metadata.json",
            "validation": "validation.json",
            "marker": "validation.ok",
            "eligible_units": "eligible_testing_units.bed.gz",
            "calibration_mask": "calibration_regions.bed.gz",
            "excluded_alleles": "excluded_alleles.tsv.gz",
            "excluded_units": "excluded_units.tsv.gz",
        }
        config["bindings"]["reference"]["lookup_source"] = source
        self.assertEqual([], validate_run_config(config))
        del source["excluded_units"]
        self.assertTrue(validate_run_config(config))

    def test_generation_rejects_indel_enabled_model(self) -> None:
        config = self.merged_manifest("tests/null_distribution/generate.run.yaml")
        config["model"]["include_indels"] = True
        self.assertTrue(validate_run_config(config))

    def test_unknown_workflow_mode_is_rejected(self) -> None:
        config = self.merged_manifest("tests/null_distribution/generate.run.yaml")
        config["workflow"]["mode"] = "unknown"
        self.assertTrue(validate_run_config(config))

    def test_mutation_policy_requires_a_rate_binding(self) -> None:
        config = self.merged_manifest("tests/null_distribution/mutation.run.yaml")
        self.assertEqual([], validate_run_config(config))
        del config["generation"]["mutation_rate"]
        self.assertTrue(validate_run_config(config))


if __name__ == "__main__":
    unittest.main()
