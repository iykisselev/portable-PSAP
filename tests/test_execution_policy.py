"""Contracts for named, scheduler-independent workflow execution policies."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from execution_policy import (  # noqa: E402
    EXECUTION_POLICIES,
    execution_policy,
    resource_mem_mb,
    resource_runtime_minutes,
)


class ExecutionPolicyTests(unittest.TestCase):
    def test_named_cadd_policy_preserves_the_recorded_scaling_coefficients(self) -> None:
        with TemporaryDirectory() as temporary:
            resource = Path(temporary) / "resource.bin"
            resource.write_bytes(b"x" * 1024)

            self.assertEqual(2304, resource_mem_mb("stage_cadd_snv", [resource], 1))
            self.assertEqual(
                380,
                resource_runtime_minutes("stage_cadd_snv", [resource], 1),
            )

    def test_unknown_staging_role_uses_the_explicit_default_policy(self) -> None:
        missing = Path("does-not-exist")

        self.assertEqual(4096, resource_mem_mb("stage_par_regions", [missing], 1))
        self.assertEqual(
            360,
            resource_runtime_minutes("stage_par_regions", [missing], 2),
        )
        with self.assertRaisesRegex(ValueError, "Unknown execution policy"):
            execution_policy("unknown")

    def test_policy_registry_is_immutable(self) -> None:
        with self.assertRaises(TypeError):
            EXECUTION_POLICIES["new"] = execution_policy("stage_vcf")  # type: ignore[index]


if __name__ == "__main__":
    unittest.main()
