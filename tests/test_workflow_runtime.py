"""Unit tests for scheduler-independent resource scaling."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from workflow_runtime import input_size_mb, scaled_mem_mb, scaled_runtime_minutes  # noqa: E402


class WorkflowRuntimeTests(unittest.TestCase):
    def test_missing_inputs_do_not_change_baseline_or_retry_scaling(self) -> None:
        missing = Path("does-not-exist")

        self.assertEqual(0.0, input_size_mb([missing]))
        self.assertEqual(4096, scaled_mem_mb([missing], 1, 4096, 512))
        self.assertEqual(240, scaled_runtime_minutes([missing], 2, 120, 30))

    def test_declared_input_size_contributes_whole_gibibyte_increments(self) -> None:
        with TemporaryDirectory() as temporary:
            resource = Path(temporary) / "resource.bin"
            resource.write_bytes(b"x" * 1024)

            self.assertGreater(input_size_mb([resource]), 0)
            self.assertEqual(4608, scaled_mem_mb([resource], 1, 4096, 512))
            self.assertEqual(150, scaled_runtime_minutes([resource], 1, 120, 30))


if __name__ == "__main__":
    unittest.main()
