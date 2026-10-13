"""Columnar SNV membership construction coverage."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import polars as pl
from polars.testing import assert_frame_equal


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_partitions import GenomePartition  # noqa: E402
from generation_reference import ReferenceReader, membership_rows  # noqa: E402


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
    def test_nested_and_unsorted_masks_do_not_omit_callable_bases(self):
        for mask in ([(0, 10), (3, 5)], [(3, 5), (0, 10)]):
            with self.subTest(mask=mask):
                rows = membership_rows(
                    GenomePartition("p", "1", 0, 10, 0),
                    mask,
                    [("unit", 4, 10)],
                    StaticReference("A" * 10),
                )
                self.assertEqual(
                    list(range(5, 11)), rows["POS"].unique().sort().to_list()
                )
                self.assertEqual(18, rows.height)

    def test_reference_reader_uses_custom_fai_without_creating_a_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta, fai = Path(tmp) / "ref.fa", Path(tmp) / "custom.fai"
            fasta.write_text(">chr1\nACGT\n")
            fai.write_text("chr1\t4\t6\t4\t5\n")
            reader = ReferenceReader(fasta, fai)
            try:
                self.assertEqual("ACGT", reader.fetch("1", 0, 4))
                self.assertFalse(Path(str(fasta) + ".fai").exists())
            finally:
                reader.close()

    def test_fallback_reader_rejects_truncated_fasta(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            patch("generation_reference.pysam", None),
        ):
            fasta, fai = Path(tmp) / "ref.fa", Path(tmp) / "ref.fai"
            fasta.write_text(">chr1\nAC")
            fai.write_text("chr1\t4\t6\t4\t5\n")
            reader = ReferenceReader(fasta, fai)
            with self.assertRaisesRegex(ValueError, "truncated"):
                reader.fetch("1", 0, 4)

    def test_memberships_are_columnar_but_preserve_base_and_alternate_order(
        self,
    ) -> None:
        rows = membership_rows(
            GenomePartition("1_000000000000_000000000005", "1", 0, 5, 0),
            mask=[(0, 5)],
            intervals=[("unit", 0, 5)],
            reference=StaticReference("ACGTN"),
        )

        expected = pl.DataFrame(
            [
                ("unit", "1", 1, "A", "C"),
                ("unit", "1", 1, "A", "G"),
                ("unit", "1", 1, "A", "T"),
                ("unit", "1", 2, "C", "A"),
                ("unit", "1", 2, "C", "G"),
                ("unit", "1", 2, "C", "T"),
                ("unit", "1", 3, "G", "A"),
                ("unit", "1", 3, "G", "C"),
                ("unit", "1", 3, "G", "T"),
                ("unit", "1", 4, "T", "A"),
                ("unit", "1", 4, "T", "C"),
                ("unit", "1", 4, "T", "G"),
            ],
            schema={
                "testing_unit_id": pl.String,
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
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
