"""Contracts for named, scheduler-independent workflow execution policies."""

from __future__ import annotations

from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from execution_policy import (  # noqa: E402
    EXECUTION_POLICIES,
    RULE_POLICY_NAMES,
    TABLE_RESOURCE_POLICY_NAMES,
    execution_policy,
    policy_resources,
    resource_mem_mb,
    resource_runtime_minutes,
    table_resource_policy,
)


class ExecutionPolicyTests(unittest.TestCase):
    def test_declared_policies_cover_every_fixed_rule_and_staging_role(self) -> None:
        expected_rule_names = {
            "prepare_generation_scope",
            "prepare_generation_partition",
            "calculate_partition_null_distribution",
            "reduce_null_distribution_bucket",
            "finalize_null_distribution_bundle",
            "validate_null_distribution_bundle",
            "validate_cohort",
            "validate_resources",
            "validate_model",
            "merge_validation",
            "require_valid_inputs",
            "extract_gt",
            "standardize_gt",
            "apply_sex_model",
            "join_sample_annotations",
            "annotate_variants",
            "summarize_annotations",
            "require_valid_annotations",
            "prepare_testing_units",
            "filter_calibration_regions",
            "filter_excluded_records",
            "build_sample_candidates",
            "build_control_reference",
            "build_cohort_report",
            "prepare_vcf_resource",
            "combine_exclusions",
            "apply_psap",
        }
        expected_staging_roles = {
            "calibration_regions": "stage_calibration_regions",
            "testing_units": "stage_testing_units",
            "cadd_snv": "stage_cadd_snv",
            "cadd_indel": "stage_cadd_indel",
            "variants_exclude": "stage_variants_exclude",
            "par_regions": "stage_par_regions",
        }

        self.assertEqual(RULE_POLICY_NAMES, expected_rule_names)
        self.assertEqual(dict(TABLE_RESOURCE_POLICY_NAMES), expected_staging_roles)
        self.assertEqual(
            set(EXECUTION_POLICIES),
            expected_rule_names | set(expected_staging_roles.values()),
        )

        declared_rule_names = set()
        for rule_file in (ROOT / "workflow" / "rules").rglob("*.smk"):
            for line in rule_file.read_text(encoding="utf-8").splitlines():
                match = re.match(r"(?:rule|checkpoint) ([A-Za-z_][A-Za-z0-9_]*):", line)
                if match:
                    declared_rule_names.add(match.group(1))
        self.assertEqual(
            declared_rule_names,
            expected_rule_names | {"prepare_table_resource"},
        )

    def test_policy_values_preserve_each_rule_resource_request(self) -> None:
        expected = {
            "prepare_generation_scope": (4096, 120, 0, 0, False),
            "prepare_generation_partition": (16384, 360, 0, 0, False),
            "calculate_partition_null_distribution": (8192, 240, 0, 0, False),
            "reduce_null_distribution_bucket": (4096, 120, 0, 0, False),
            "finalize_null_distribution_bundle": (2048, 120, 0, 0, False),
            "validate_null_distribution_bundle": (1024, 120, 0, 0, False),
            "validate_cohort": (4096, 30, 0, 0, False),
            "validate_resources": (4096, 30, 0, 0, False),
            "validate_model": (512, 10, 0, 0, False),
            "merge_validation": (512, 5, 0, 0, False),
            "require_valid_inputs": (512, 5, 0, 0, False),
            "extract_gt": (4096, 120, 0, 0, False),
            "standardize_gt": (8192, 120, 0, 0, False),
            "apply_sex_model": (8192, 120, 0, 0, False),
            "join_sample_annotations": (16000, 180, 0, 0, False),
            "annotate_variants": (32000, 720, 512, 30, True),
            "summarize_annotations": (1024, 30, 0, 0, False),
            "require_valid_annotations": (512, 5, 0, 0, False),
            "prepare_testing_units": (16000, 180, 0, 0, False),
            "filter_calibration_regions": (4096, 180, 0, 0, False),
            "filter_excluded_records": (4096, 180, 0, 0, False),
            "build_sample_candidates": (4096, 60, 0, 0, False),
            "build_control_reference": (16000, 360, 256, 20, True),
            "build_cohort_report": (16000, 360, 256, 20, True),
            "prepare_vcf_resource": (4096, 180, 0, 0, True),
            "combine_exclusions": (12000, 180, 0, 0, False),
            "apply_psap": (16000, 360, 0, 0, False),
            "stage_calibration_regions": (2048, 60, 0, 0, True),
            "stage_testing_units": (2048, 60, 0, 0, True),
            "stage_cadd_snv": (2048, 360, 256, 20, True),
            "stage_cadd_indel": (2048, 360, 256, 20, True),
            "stage_variants_exclude": (4096, 180, 0, 0, True),
            "stage_par_regions": (4096, 180, 0, 0, True),
        }

        observed = {
            name: (
                policy.baseline_mem_mb,
                policy.baseline_runtime_minutes,
                policy.mem_mb_per_gib,
                policy.runtime_minutes_per_gib,
                policy.scale_with_attempt,
            )
            for name, policy in EXECUTION_POLICIES.items()
        }
        self.assertEqual(observed, expected)

    def test_size_scaled_policies_preserve_recorded_retry_scaling(self) -> None:
        with TemporaryDirectory() as temporary:
            resource = Path(temporary) / "resource.bin"
            resource.write_bytes(b"x" * 1024)

            self.assertEqual(2304, resource_mem_mb("stage_cadd_snv", [resource], 1))
            self.assertEqual(
                380,
                resource_runtime_minutes("stage_cadd_snv", [resource], 1),
            )
            self.assertEqual(4608, resource_mem_mb("stage_cadd_snv", [resource], 2))
            self.assertEqual(
                760,
                resource_runtime_minutes("stage_cadd_snv", [resource], 2),
            )

    def test_fixed_policies_do_not_gain_retry_scaling(self) -> None:
        callbacks = policy_resources("apply_psap")

        self.assertEqual(16000, callbacks.mem_mb(None, [], 2))
        self.assertEqual(360, callbacks.runtime(None, [], 2))

    def test_unknown_policies_and_staging_roles_fail_explicitly(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown execution policy"):
            execution_policy("unknown")
        with self.assertRaisesRegex(ValueError, "Unknown execution policy"):
            execution_policy("stage_not_registered")
        with self.assertRaisesRegex(ValueError, "No staging policy"):
            table_resource_policy("not_a_resource_role")

        self.assertEqual(
            "stage_par_regions",
            table_resource_policy("par_regions"),
        )

    def test_policy_registry_is_immutable(self) -> None:
        with self.assertRaises(TypeError):
            EXECUTION_POLICIES["new"] = execution_policy("prepare_vcf_resource")  # type: ignore[index]


if __name__ == "__main__":
    unittest.main()
