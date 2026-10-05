"""Contracts for bounded normalized indel catalog candidates."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_indels import (  # noqa: E402
    INDEL_MEMBERSHIP_SCHEMA,
    INDEL_REJECTION_SCHEMA,
    INDEL_SCORE_SCHEMA,
    IndelCandidateError,
    indel_candidate_memberships,
    indel_score_rows,
)


class InMemoryReference:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        if chromosome != "1" or start < 0 or end > len(self.sequence):
            raise ValueError("outside test reference")
        return self.sequence[start:end]


class GenerationIndelTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.resource = self.directory / "indels.tsv"
        self.reference = InMemoryReference("CAGTAC")

    def rows(self, content: str, *, start: int = 0, end: int = 6) -> pl.DataFrame:
        self.resource.write_text(content, encoding="utf-8")
        return indel_score_rows(
            self.resource,
            None,
            "1",
            start,
            end,
            self.reference,
            max_reference_length=3,
            max_alternate_length=3,
        )

    def test_reader_uses_anchor_owner_window_and_preserves_duplicates(self) -> None:
        rows = self.rows(
            "chr1\t2\tA\tAT\t0\t15\n"
            "1\t2\tA\tAT\t0\t15\n"
            "1\t7\tC\tCA\t0\t21\n"
        )
        self.assertEqual(rows.schema, INDEL_SCORE_SCHEMA)
        self.assertEqual(rows.height, 2)
        self.assertEqual(rows["POS"].to_list(), [2, 2])
        self.assertEqual(rows["SPAN_START"].to_list(), [1, 1])
        self.assertEqual(rows["SPAN_END"].to_list(), [2, 2])

    def test_reader_rejects_invalid_catalog_rows(self) -> None:
        cases = {
            "not_indel": ("1\t2\tA\tG\t0\t15\n", "not a length-changing"),
            "symbolic": ("1\t2\tA\t<DEL>\t0\t15\n", "SNV or length-changing"),
            "not_minimal": ("1\t2\tAGT\tAT\t0\t15\n", "not minimal"),
            "not_left_aligned": ("1\t2\tA\tAA\t0\t15\n", "not left-aligned"),
            "reference_mismatch": ("1\t2\tC\tCT\t0\t15\n", "does not match"),
            "length": ("1\t2\tAGTA\tA\t0\t15\n", "exceeds configured"),
            "score": ("1\t2\tA\tAT\t0\tnan\n", "non-finite"),
        }
        for name, (content, message) in cases.items():
            with self.subTest(name=name), self.assertRaisesRegex(IndelCandidateError, message):
                self.rows(content)

    def test_full_span_mask_and_unit_membership_are_distinct(self) -> None:
        candidates = pl.DataFrame(
            [
                ("1", 10, "AAA", "A", 9, 12, 20.0),
                ("1", 7, "AA", "A", 6, 8, 21.0),
                ("1", 9, "A", "AT", 8, 9, 22.0),
            ],
            schema=INDEL_SCORE_SCHEMA,
            orient="row",
        )
        memberships, rejected = indel_candidate_memberships(
            candidates,
            mask=[(8, 13)],
            unit_intervals=[("unit-a", 10, 11), ("unit-b", 11, 13)],
        )
        self.assertEqual(memberships.schema, INDEL_MEMBERSHIP_SCHEMA)
        self.assertEqual(
            [(row["testing_unit_id"], row["POS"], row["variant_class"])
             for row in memberships.sort("testing_unit_id").iter_rows(named=True)],
            [("unit-a", 10, "indel"), ("unit-b", 10, "indel")],
        )
        self.assertEqual(rejected.schema, INDEL_REJECTION_SCHEMA)
        self.assertEqual(
            {(row["POS"], row["rejection_reason"]) for row in rejected.iter_rows(named=True)},
            {
                (7, "span_not_contained_in_calibration_mask"),
                (9, "no_eligible_testing_unit_overlap"),
            },
        )


if __name__ == "__main__":
    unittest.main()
