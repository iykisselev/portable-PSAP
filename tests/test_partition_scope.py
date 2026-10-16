"""Contracts for biological partition scope accounting."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partition_scope import build_partition_scope  # noqa: E402


class PartitionScopeTests(unittest.TestCase):
    def test_precedence_is_applied_before_retention(self) -> None:
        scope = build_partition_scope(
            [
                ("cds", {"GENE": [("1", 100, 200)]}),
                ("intronic", {"GENE": [("1", 50, 250)]}),
            ],
            {"1": [(0, 1000)]},
            {"1": [(50, 100), (200, 225)]},
            0.5,
        )
        cds, intronic = scope.rows
        self.assertEqual("fully_excluded", cds["status"])
        self.assertEqual(100, cds["raw_bases"])
        self.assertEqual("eligible", intronic["status"])
        self.assertEqual(200, intronic["raw_bases"])
        self.assertEqual(100, intronic["priority_removed_bases"])
        self.assertEqual(100, intronic["resolved_bases"])
        self.assertEqual(75, intronic["retained_bases"])

    def test_empty_after_priority_is_distinct_from_mask_exclusion(self) -> None:
        scope = build_partition_scope(
            [("first", {"G": [("1", 0, 10)]}), ("second", {"G": [("1", 0, 10)]})],
            {"1": [(0, 10)]}, {"1": [(0, 10)]}, 0.5,
        )
        self.assertEqual("empty_after_priority", scope.rows[1]["status"])


if __name__ == "__main__":
    unittest.main()
