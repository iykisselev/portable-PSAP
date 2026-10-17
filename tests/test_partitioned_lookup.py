"""Partition/model public lookup reduction contracts."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partitioned_lookup import partition_model_lookup_rows  # noqa: E402


class PartitionedLookupTests(unittest.TestCase):
    def test_every_model_has_a_separate_partition_table(self):
        statistics = pl.DataFrame(
            {
                "testing_unit_id": ["AUTO", "XGENE"],
                "partition_id": ["cds", "cds"],
                "role": ["het", "het"],
                "bin": [1, 1],
                "certain_count": [0, 0],
                "log_none_sum": [-0.1, -0.2],
                "event_count": [1, 1],
                "snv_event_count": [1, 1],
                "indel_event_count": [0, 0],
            }
        )
        tables, _ = partition_model_lookup_rows(
            statistics,
            [
                {"testing_unit_id": "AUTO", "partition_id": "cds", "CHROM": "1"},
                {"testing_unit_id": "XGENE", "partition_id": "cds", "CHROM": "X"},
            ],
            hemizygote=True,
            bin_count=2,
        )
        self.assertEqual(5, len(tables))
        self.assertEqual(["AUTO"], [row[0] for row in tables[("cds", "DOM-het")]])
        self.assertEqual(["XGENE"], [row[0] for row in tables[("cds", "X-linked-het")]])
        self.assertEqual([], tables[("cds", "REC-hom")])


if __name__ == "__main__":
    unittest.main()
