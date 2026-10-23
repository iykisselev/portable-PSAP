"""Contracts for bounded, normalized indel candidate generation."""

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
    partitioned_indel_candidate_memberships,
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
        self.resource = Path(temporary.name) / "indels.tsv"
        self.reference = InMemoryReference("CAGTAC")

    def test_reader_uses_owner_window_and_preserves_duplicate_catalog_rows(
        self,
    ) -> None:
        self.resource.write_text(
            "chr1\t2\tA\tAT\t0\t15\n1\t2\tA\tAT\t0\t15\n1\t7\tC\tCA\t0\t21\n",
            encoding="utf-8",
        )
        rows = indel_score_rows(
            self.resource,
            None,
            "1",
            0,
            6,
            self.reference,
            max_reference_length=3,
            max_alternate_length=3,
        )
        self.assertEqual(INDEL_SCORE_SCHEMA, rows.schema)
        self.assertEqual([2, 2], rows["POS"].to_list())

    def test_reader_rejects_invalid_catalog_rows(self) -> None:
        self.resource.write_text("1\t2\tA\tG\t0\t15\n", encoding="utf-8")
        with self.assertRaisesRegex(IndelCandidateError, "not a length-changing"):
            indel_score_rows(
                self.resource,
                None,
                "1",
                0,
                6,
                self.reference,
                max_reference_length=3,
                max_alternate_length=3,
            )

    def test_membership_uses_full_ref_span_and_unit_overlap(self) -> None:
        candidates = pl.DataFrame(
            [("1", 10, "AAA", "A", 9, 12, 20.0), ("1", 7, "AA", "A", 6, 8, 21.0)],
            schema=INDEL_SCORE_SCHEMA,
            orient="row",
        )
        memberships, rejected = indel_candidate_memberships(
            candidates,
            mask=[(8, 13)],
            unit_intervals=[("unit-a", 10, 11), ("unit-b", 11, 13)],
        )
        self.assertEqual(INDEL_MEMBERSHIP_SCHEMA, memberships.schema)
        self.assertEqual({"unit-a", "unit-b"}, set(memberships["testing_unit_id"]))
        self.assertEqual(INDEL_REJECTION_SCHEMA, rejected.schema)
        self.assertEqual(
            ["span_not_contained_in_calibration_mask"],
            rejected["rejection_reason"].to_list(),
        )

    def test_partitioned_membership_uses_the_earliest_overlapping_owner(self) -> None:
        candidates = pl.DataFrame(
            [("1", 10, "AAA", "A", 9, 12, 20.0)],
            schema=INDEL_SCORE_SCHEMA,
            orient="row",
        )
        memberships, _ = partitioned_indel_candidate_memberships(
            candidates,
            mask=[(8, 13)],
            intervals=[
                ("unit-a", "cds", 0, 10, 11),
                ("unit-a", "intronic", 1, 11, 13),
            ],
        )
        self.assertEqual(["cds"], memberships["partition_id"].to_list())


if __name__ == "__main__":
    unittest.main()
