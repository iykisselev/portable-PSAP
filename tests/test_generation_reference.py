"""Columnar SNV membership construction coverage."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import polars as pl
from polars.testing import assert_frame_equal


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_partitions import GenomePartition  # noqa: E402
from generation_reference import membership_rows  # noqa: E402


class StaticReference:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        self.assertEqual(chromosome, "1")
        return self.sequence[start:end]

    def assertEqual(self, left, right) -> None:
        if left != right:
            raise AssertionError(f"{left!r} != {right!r}")


class GenerationReferenceTests(unittest.TestCase):
    def test_memberships_are_columnar_but_preserve_base_and_alternate_order(self) -> None:
        rows = membership_rows(
            GenomePartition("1_000000000000_000000000005", "1", 0, 5, 0),
            mask=[(0, 5)],
            intervals=[("unit", 0, 5)],
            reference=StaticReference("ACGTN"),
        )

        expected = pl.DataFrame(
            [
                ("unit", "1", 1, "A", "C"), ("unit", "1", 1, "A", "G"), ("unit", "1", 1, "A", "T"),
                ("unit", "1", 2, "C", "A"), ("unit", "1", 2, "C", "G"), ("unit", "1", 2, "C", "T"),
                ("unit", "1", 3, "G", "A"), ("unit", "1", 3, "G", "C"), ("unit", "1", 3, "G", "T"),
                ("unit", "1", 4, "T", "A"), ("unit", "1", 4, "T", "C"), ("unit", "1", 4, "T", "G"),
            ],
            schema={
                "testing_unit_id": pl.String, "CHROM": pl.String, "POS": pl.Int64,
                "REF": pl.String, "ALT": pl.String,
            },
            orient="row",
        )
        assert_frame_equal(rows, expected)

    def test_memberships_remain_distinct_for_overlapping_direct_inputs(self) -> None:
        rows = membership_rows(
            GenomePartition("1_000000000000_000000000002", "1", 0, 2, 0),
            mask=[(0, 2)],
            intervals=[("unit", 0, 2), ("unit", 1, 2)],
            reference=StaticReference("AA"),
        )

        self.assertEqual(rows.height, 6)


if __name__ == "__main__":
    unittest.main()
