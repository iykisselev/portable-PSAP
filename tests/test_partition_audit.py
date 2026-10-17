"""Testing-unit ID-overlap audit contracts."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partition_audit import partition_id_overlap_rows  # noqa: E402


class PartitionAuditTests(unittest.TestCase):
    def test_pair_metrics_and_low_overlap_warning(self) -> None:
        row = partition_id_overlap_rows(
            [("cds", {"A": None, "B": None}), ("intronic", {"B": None, "C": None})]
        )[0]
        self.assertEqual(1, row["intersection_count"])
        self.assertEqual(3, row["union_count"])
        self.assertEqual(0.5, row["minimum_overlap_fraction"])
        self.assertFalse(row["low_overlap_warning"])

    def test_empty_sets_are_explicit(self) -> None:
        row = partition_id_overlap_rows([("cds", {}), ("intronic", {"A": None})])[0]
        self.assertEqual("left_empty", row["empty_set_status"])
        self.assertIsNone(row["minimum_overlap_fraction"])


if __name__ == "__main__":
    unittest.main()
