"""Unit exclusion for immutable workflow mode and binding resolution."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import sys
import tempfile
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
from bundle_contract import INDEL_GENERATED_BUNDLE_CONTRACT  # noqa: E402


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

    def test_nonfinite_score_max_is_rejected_before_scoring(self):
        for value in (float("inf"), float("nan")):
            with self.subTest(value=value):
                manifest = self.manifest("tests/smoke/run.yaml")
                manifest["model"]["score_max"] = value
                with self.assertRaisesRegex(
                    WorkflowPlanError, "model.score_max.*finite"
                ):
                    build_workflow_plan(manifest)

    def directory_manifest(
        self, directory, *, contract="generated_bundle_v2", hem=True
    ):
        """Small declarations stand in for component bytes at the planning boundary."""
        from bundle_contract import BUNDLE_COMPONENT_FILENAMES

        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["model"].update(hemizygote=hem, include_indels=contract.endswith("v3"))
        reference = manifest["bindings"]["reference"]
        for name in ("testing_unit_bed", "exclusion", "exclusions"):
            reference.pop(name, None)
        reference["lookup_source"] = {"type": "bundle", "directory": str(directory)}
        roles = ["het", "hom"] + (["hem_chrX"] if hem else [])
        components = {
            name: {
                "path": BUNDLE_COMPONENT_FILENAMES[
                    "eligible_bed" if name == "eligible_testing_units" else name
                ]
            }
            for name in (
                "eligible_testing_units",
                "calibration_regions",
                "excluded_alleles",
                "excluded_units",
                *(f"lookup_{role}" for role in roles),
            )
        }
        for component in components.values():
            (directory / component["path"]).touch()
        for name in ("validation", "marker"):
            (directory / BUNDLE_COMPONENT_FILENAMES[name]).touch()
        metadata = {
            "contract": contract,
            "binning": {"score_max": 23, "bin_width": 0.2},
            "model_assumptions": {"roles": roles},
            "components": components,
        }
        (directory / "metadata.json").write_text(json.dumps(metadata))
        return manifest, metadata

    def test_directory_resolves_declared_grid_and_scope_for_all_reader_versions(self):
        for contract in (
            "generated_bundle_v1",
            "generated_bundle_v2",
            "generated_bundle_v3",
        ):
            for hem in (False, True):
                with (
                    self.subTest(contract=contract, hem=hem),
                    tempfile.TemporaryDirectory() as tmp,
                ):
                    directory = Path(tmp)
                    manifest, _ = self.directory_manifest(
                        directory, contract=contract, hem=hem
                    )
                    original = copy.deepcopy(manifest)
                    plan = build_workflow_plan(manifest)
                    lookup = plan.apply.lookup_source
                    self.assertEqual(23, lookup.score_max)
                    self.assertEqual(0.2, lookup.bin_width)
                    self.assertEqual(contract, lookup.contract)
                    self.assertEqual(hem, "hem_chrX" in lookup.tables)
                    self.assertEqual(
                        str(directory / "eligible_testing_units.bed.gz"),
                        plan.apply.effective_reference.testing_unit_bed.data,
                    )
                    self.assertIsNone(
                        plan.apply.configured_reference.testing_unit_bed.data
                    )
                    self.assertEqual(original, manifest)
                    if contract == "generated_bundle_v3" and hem:
                        example = self.manifest("runs/apply.bundle.example.yaml")
                        example["bindings"]["reference"]["lookup_source"][
                            "directory"
                        ] = str(directory)
                        self.assertEqual(
                            lookup, build_workflow_plan(example).apply.lookup_source
                        )

    def test_directory_and_explicit_forms_resolve_equivalently(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest, _ = self.directory_manifest(directory)
            directory_plan = build_workflow_plan(manifest)
            lookup = directory_plan.apply.lookup_source
            manifest["model"]["score_max"] = lookup.score_max
            manifest["generation_settings"]["score_bin_width"] = lookup.bin_width
            manifest["bindings"]["reference"]["lookup_source"] = {
                "type": "bundle",
                "contract": lookup.contract,
                "tables": dict(lookup.tables),
                **{
                    name: getattr(lookup, name)
                    for name in (
                        "metadata",
                        "validation",
                        "marker",
                        "eligible_units",
                        "calibration_mask",
                        "excluded_units",
                    )
                },
                "excluded_alleles": lookup.excluded_alleles[0],
            }
            self.assertEqual(lookup, build_workflow_plan(manifest).apply.lookup_source)

    def test_directory_rejects_mixed_forms_and_reports_selected_branch(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, _ = self.directory_manifest(Path(tmp))
            source = manifest["bindings"]["reference"]["lookup_source"]
            source["contract"] = "generated_bundle_v2"
            with self.assertRaisesRegex(WorkflowPlanError, "do not combine") as caught:
                build_workflow_plan(manifest)
            self.assertNotIn("external_legacy_v0", str(caught.exception))

    def test_directory_rejects_incomplete_invalid_or_escaping_declarations(self):
        for failure in (
            "metadata",
            "component",
            "certificate",
            "contract",
            "grid",
            "escape",
        ):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                manifest, metadata = self.directory_manifest(directory)
                if failure == "metadata":
                    (directory / "metadata.json").unlink()
                elif failure == "component":
                    (directory / "lookup_het.tsv.gz").unlink()
                elif failure == "certificate":
                    (directory / "validation.ok").unlink()
                else:
                    if failure == "contract":
                        metadata["contract"] = "unknown"
                    elif failure == "grid":
                        metadata["binning"]["bin_width"] = 0
                    else:
                        metadata["components"]["lookup_het"]["path"] = (
                            "../lookup_het.tsv.gz"
                        )
                    (directory / "metadata.json").write_text(json.dumps(metadata))
                with self.assertRaisesRegex(
                    WorkflowPlanError, "lookup_source.directory"
                ):
                    build_workflow_plan(manifest)

    def test_directory_rejects_hem_mismatch_before_rule_input_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, _ = self.directory_manifest(Path(tmp), hem=False)
            manifest["model"]["hemizygote"] = True
            with self.assertRaisesRegex(WorkflowPlanError, "model.hemizygote"):
                build_workflow_plan(manifest)

    def test_directory_supports_invocation_relative_paths_and_symlinked_components(
        self,
    ):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest, _ = self.directory_manifest(directory)
            component = directory / "lookup_het.tsv.gz"
            component.unlink()
            target = directory / "source.tsv.gz"
            target.touch()
            component.symlink_to(target)
            relative = os.path.relpath(directory)
            manifest["bindings"]["reference"]["lookup_source"]["directory"] = relative
            plan = build_workflow_plan(manifest)
            self.assertEqual(
                str(Path(relative) / component.name),
                plan.apply.lookup_source.tables["het"],
            )

    def test_bundle_fields_outside_lookup_source_have_a_migration_hint(self):
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["bindings"]["reference"]["metadata"] = "bundle/metadata.json"
        with self.assertRaisesRegex(
            WorkflowPlanError, "inside bindings.reference.lookup_source"
        ):
            build_workflow_plan(manifest)

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
            "tests/null_distribution/exclusion.bed",
            plan.apply.configured_reference.exclusion.data,
        )
        self.assertEqual(
            str(plan.generation.bundle_dir / "calibration_regions.bed.gz"),
            plan.apply.effective_reference.calibration_regions.data,
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

    def test_generation_plan_exposes_staged_indel_contract_bindings(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["generation"]["score_indel"] = {
            "data": "tests/smoke/data/cadd_indel.tsv",
            "index": None,
        }
        manifest["generation"]["mutation_rate_indel"] = {
            "data": "tests/null_distribution/mutation_rates.tsv",
            "index": None,
        }
        plan = build_workflow_plan(manifest)

        assert plan.generation is not None
        self.assertEqual(
            "tests/smoke/data/cadd_indel.tsv", plan.generation.score_indel.data
        )
        self.assertEqual(
            "tests/null_distribution/mutation_rates.tsv",
            plan.generation.mutation_rate_indel.data,
        )
        settings = plan.generation_settings
        self.assertEqual("score_resource", settings.indel_candidate_universe)
        self.assertEqual("reference_span_overlap", settings.indel_membership)
        self.assertEqual("reference_span_contained", settings.indel_mask_policy)
        self.assertEqual(50, settings.indel_max_reference_length)
        self.assertEqual(50, settings.indel_max_alternate_length)

    def test_generation_rejects_indel_variant_class_before_execution_is_ready(
        self,
    ) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["generation_settings"]["variant_classes"] = ["snv", "indel"]
        with self.assertRaisesRegex(WorkflowPlanError, "variant_classes must be"):
            build_workflow_plan(manifest)

    def test_generate_plan_activates_the_combined_indel_contract(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["model"]["include_indels"] = True
        manifest["generation_settings"]["variant_classes"] = ["snv", "indel"]
        manifest["generation"]["score_indel"] = {
            "data": "tests/smoke/data/cadd_indel.tsv",
            "index": None,
        }
        plan = build_workflow_plan(manifest)

        assert plan.generation is not None
        self.assertEqual(
            "tests/smoke/data/cadd_indel.tsv", plan.generation.score_indel.data
        )

    def test_full_plan_uses_v3_lookup_contract_for_combined_generation(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["model"]["include_indels"] = True
        manifest["generation_settings"]["variant_classes"] = ["snv", "indel"]
        manifest["bindings"]["reference"]["cadd"]["indel"] = {
            "data": "tests/smoke/data/cadd_indel.tsv",
            "index": None,
        }
        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        self.assertEqual(
            INDEL_GENERATED_BUNDLE_CONTRACT, plan.apply.lookup_source.contract
        )

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

    def test_optional_apply_exclusion_remains_absent(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["bindings"]["reference"]["exclusion"] = {
            "data": None,
            "index": None,
        }
        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        self.assertIsNone(plan.apply.effective_reference.exclusion.data)
        self.assertIsNone(plan.apply.effective_reference.calibration_regions)
        self.assertIsNone(plan.apply.lookup_source.calibration_mask)

    def test_production_generation_requires_each_partition_resource_index(self) -> None:
        for mode in ("generate", "full"):
            for missing in ("population", "scores", "mutation"):
                with self.subTest(mode=mode, missing=missing):
                    manifest = self.manifest(f"tests/null_distribution/{mode}.run.yaml")
                    manifest["generation"]["execution"].pop(
                        "allow_unindexed_fixture_inputs"
                    )
                    population = manifest["generation"]["population_vcf"]
                    scores = (
                        manifest["generation"]["score_snv"]
                        if mode == "generate"
                        else manifest["bindings"]["reference"]["cadd"]["snv"]
                    )
                    mutation = {
                        "data": "/production/rates.tsv.gz",
                        "index": "/production/rates.tsv.gz.csi",
                    }
                    manifest["generation"]["mutation_rate"] = mutation
                    manifest["generation"]["unobserved_policy"] = "mutation_rate"
                    population["index"] = "/production/population.vcf.gz.tbi"
                    scores["index"] = "/production/scores.tsv.gz.csi"
                    resources = {
                        "population": population,
                        "scores": scores,
                        "mutation": mutation,
                    }
                    # Indexed resources need not exist at plan construction time.
                    build_workflow_plan(manifest)
                    resources[missing]["index"] = None
                    key = {
                        "population": "generation.population_vcf",
                        "mutation": "generation.mutation_rate",
                        "scores": "generation.score_snv"
                        if mode == "generate"
                        else "bindings.reference.cadd.snv",
                    }[missing]
                    with self.assertRaises(WorkflowPlanError) as raised:
                        build_workflow_plan(manifest)
                    self.assertIn(f"{key}.index is required", str(raised.exception))
                    self.assertIn(
                        "each partition scans the entire file", str(raised.exception)
                    )

    def test_fixture_exception_never_allows_unindexed_compressed_inputs(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["generation"]["population_vcf"]["data"] = "/production/gnomad.vcf.gz"
        with self.assertRaisesRegex(
            WorkflowPlanError, "generation.population_vcf.index is required"
        ):
            build_workflow_plan(manifest)

    def test_fixture_exception_requires_boolean_not_truthy_text(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["generation"]["execution"]["allow_unindexed_fixture_inputs"] = "false"
        with self.assertRaisesRegex(WorkflowPlanError, "must be a boolean"):
            build_workflow_plan(manifest)

    def test_model_flags_require_real_booleans(self) -> None:
        for field in ("hemizygote", "include_indels"):
            for invalid_value in ("false", "true", 0, 1, None):
                with self.subTest(field=field, invalid_value=invalid_value):
                    manifest = self.manifest("tests/smoke/run.yaml")
                    manifest["model"][field] = invalid_value

                    with self.assertRaisesRegex(
                        WorkflowPlanError, f"model.{field} must be a boolean"
                    ):
                        build_workflow_plan(manifest)

    def test_full_plan_requires_one_shared_par_role(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["reference"]["exclusions"]["par-regions"]["role"] = "variants_exclude"

        with self.assertRaisesRegex(WorkflowPlanError, "exactly one.*par_regions"):
            build_workflow_plan(manifest)

    def test_exclusions_are_resolved_by_identifier_not_manifest_order(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        bindings = manifest["bindings"]["reference"]["exclusions"]
        manifest["bindings"]["reference"]["exclusions"] = {
            "par-regions": bindings["par-regions"],
            "encode-hg38-blacklist-v2": bindings["encode-hg38-blacklist-v2"],
        }

        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        resolved = {
            item.identifier: item.binding.data
            for item in plan.apply.effective_exclusions
        }
        self.assertEqual(
            {
                "encode-hg38-blacklist-v2": "tests/smoke/data/exclude_full_span.bed",
                "par-regions": "tests/smoke/data/par_regions.bed",
            },
            resolved,
        )

    def test_exclusion_identifier_mismatches_fail_before_resource_registry(
        self,
    ) -> None:
        missing = self.manifest("tests/smoke/run.yaml")
        del missing["bindings"]["reference"]["exclusions"]["encode-hg38-blacklist-v2"]
        with self.assertRaisesRegex(
            WorkflowPlanError, "missing active bindings.*encode-hg38-blacklist-v2"
        ):
            build_workflow_plan(missing)

        unknown = self.manifest("tests/smoke/run.yaml")
        unknown["bindings"]["reference"]["exclusions"]["unexpected"] = {
            "data": "unexpected.bed",
            "index": None,
        }
        with self.assertRaisesRegex(
            WorkflowPlanError, "no scientific metadata.*unexpected"
        ):
            build_workflow_plan(unknown)

    def test_disabled_scientific_exclusion_does_not_require_or_stage_a_binding(
        self,
    ) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["reference"]["exclusions"]["encode-hg38-blacklist-v2"]["enabled"] = (
            False
        )

        plan = build_workflow_plan(manifest)

        assert plan.apply is not None
        self.assertEqual(
            ("par-regions",),
            tuple(item.identifier for item in plan.apply.effective_exclusions),
        )

    def test_generation_execution_defaults_remain_available_to_direct_callers(
        self,
    ) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["generation"]["execution"] = {"allow_unindexed_fixture_inputs": True}

        plan = build_workflow_plan(manifest)

        assert plan.generation is not None
        self.assertEqual(5_000_000, plan.generation.partition_size_bp)
        self.assertEqual(250, plan.generation.unit_bucket_size)
        self.assertEqual("zstd", plan.generation.parquet_compression)

    def test_direct_plan_validation_rejects_unknown_configuration_fields(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["bindings"]["reference"]["mistyped_resource"] = {}

        with self.assertRaisesRegex(
            WorkflowPlanError, "Additional properties are not allowed"
        ):
            build_workflow_plan(manifest)

    def test_legacy_exclusion_lists_explain_the_migration(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["bindings"]["reference"]["exclusions"] = [
            {"data": "old-style.bed", "index": None}
        ]

        with self.assertRaisesRegex(WorkflowPlanError, "mapping keyed by exclusion ID"):
            build_workflow_plan(manifest)

    def test_legacy_normalization_warns_without_changing_the_source_mapping(
        self,
    ) -> None:
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

    def assert_rejected(self, manifest: dict) -> None:
        with self.assertRaises(WorkflowPlanError):
            build_workflow_plan(manifest)

    def test_science_configuration_rejects_removed_reference_metadata(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["reference"]["lookup_tables"] = {}
        self.assert_rejected(manifest)

        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["reference"]["exclusion"]["optional"] = True
        self.assert_rejected(manifest)

    def test_unknown_workflow_mode_is_rejected(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["workflow"]["mode"] = "unknown"
        self.assert_rejected(manifest)

    def test_generate_rejects_apply_inputs(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["inputs"] = {
            "vcf": {"data": "cohort.vcf", "index": None},
            "ped": "cohort.ped",
            "samples": None,
        }
        self.assert_rejected(manifest)

    def test_generation_rejects_indel_enabled_model(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        manifest["model"]["include_indels"] = True
        self.assert_rejected(manifest)

    def test_sex_specific_hem_generation_requires_male_af_tag(self) -> None:
        manifest = self.manifest("tests/null_distribution/generate.run.yaml")
        del manifest["generation"]["af_source"]["male"]
        self.assert_rejected(manifest)

    def test_mutation_policy_requires_a_rate_binding(self) -> None:
        manifest = self.manifest("tests/null_distribution/mutation.run.yaml")
        del manifest["generation"]["mutation_rate"]
        self.assert_rejected(manifest)

    def test_full_rejects_a_second_shared_resource_binding_set(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["generation"]["score_snv"] = manifest["bindings"]["reference"]["cadd"][
            "snv"
        ]
        self.assert_rejected(manifest)

    def test_full_rejects_external_lookup_source(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["bindings"]["reference"]["lookup_source"] = {
            "type": "external",
            "contract": "external_legacy_v0",
            "tables": {"het": "het.tsv", "hom": "hom.tsv"},
        }
        self.assert_rejected(manifest)

    def test_full_allows_no_exclusion_in_the_shared_reference_bindings(self) -> None:
        manifest = self.manifest("tests/null_distribution/full.run.yaml")
        manifest["bindings"]["reference"]["exclusion"] = {"data": None, "index": None}
        build_workflow_plan(manifest)

    def test_apply_without_hemizygote_does_not_require_hem_lookup(self) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["model"]["hemizygote"] = False
        del manifest["bindings"]["reference"]["lookup_source"]["tables"]["hem_chrX"]
        build_workflow_plan(manifest)

    def test_bundle_lookup_source_requires_every_certified_scope_component(
        self,
    ) -> None:
        manifest = self.manifest("tests/smoke/run.yaml")
        manifest["model"]["hemizygote"] = False
        manifest["model"]["include_indels"] = False
        source = {
            "type": "bundle",
            "contract": GENERATED_BUNDLE_CONTRACT,
            "tables": {"het": "lookup_het.tsv.gz", "hom": "lookup_hom.tsv.gz"},
            "metadata": "metadata.json",
            "validation": "validation.json",
            "marker": "validation.ok",
            "eligible_units": "eligible_testing_units.bed.gz",
            "calibration_mask": "calibration_regions.bed.gz",
            "excluded_alleles": "excluded_alleles.tsv.gz",
            "excluded_units": "excluded_units.tsv.gz",
        }
        manifest["bindings"]["reference"]["lookup_source"] = source
        build_workflow_plan(manifest)
        del source["excluded_units"]
        self.assert_rejected(manifest)


if __name__ == "__main__":
    unittest.main()
