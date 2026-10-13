"""Shared text, BED, contig, and FAI parsing contracts."""

from __future__ import annotations

import gzip
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from genomic_contracts import read_fai  # noqa: E402
from io_contracts import BedFormatError, iter_bed, open_text  # noqa: E402


class TextAndFaiContractTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def test_text_reader_uses_gzip_magic_not_filename_suffix(self) -> None:
        compressed = self.directory / "resource.txt"
        with gzip.open(compressed, "wt", encoding="utf-8") as handle:
            handle.write("compressed\n")
        plain = self.directory / "resource.gz"
        plain.write_text("plain\n", encoding="utf-8")

        with open_text(compressed) as handle:
            self.assertEqual("compressed\n", handle.read())
        with open_text(plain) as handle:
            self.assertEqual("plain\n", handle.read())
        with self.assertRaisesRegex(ValueError, "gzip-compressed"):
            open_text(plain, require_gzip=True)

    def test_bed_reader_skips_directives_and_normalizes_validation(self) -> None:
        bed = self.directory / "regions.bed"
        bed.write_text(
            "# comment\ntrack name=regions\nbrowser position chr1\nchr1\t2\t5\tunit\n",
            encoding="utf-8",
        )
        self.assertEqual(
            [(4, "chr1", 2, 5, "unit")],
            [
                (row.line_number, row.chromosome, row.start, row.end, row.fields[3])
                for row in iter_bed(bed, 4)
            ],
        )

    def test_bed_contig_prefixes_are_not_mistaken_for_directives(self):
        bed = self.directory / "prefixes.bed"
        bed.write_text("track1\t0\t2\nbrowser_alt\t2\t4\n")
        self.assertEqual(
            ["track1", "browser_alt"], [row.chromosome for row in iter_bed(bed)]
        )

    def test_bed_reader_can_collect_invalid_rows_or_fail_strictly(self) -> None:
        bed = self.directory / "invalid.bed"
        bed.write_text("1\t4\t2\n1\t5\t7\n", encoding="utf-8")
        errors: list[str] = []
        self.assertEqual(
            [(5, 7)],
            [
                (row.start, row.end)
                for row in iter_bed(bed, on_error=lambda exc: errors.append(str(exc)))
            ],
        )
        self.assertEqual(
            ["BED line 1 must satisfy 0 <= start < end; received 4, 2"], errors
        )
        with self.assertRaises(BedFormatError):
            list(iter_bed(bed))

    def test_fai_reader_preserves_layout_under_canonical_contig_keys(self) -> None:
        fai = self.directory / "reference.fai"
        fai.write_text("chr1\t10\t4\t5\t6\nchrM\t3\t20\t3\t4\n", encoding="utf-8")
        records = read_fai(fai)
        self.assertEqual("chr1", records["1"].name)
        self.assertEqual(
            (10, 4, 5, 6),
            (
                records["1"].length,
                records["1"].offset,
                records["1"].line_bases,
                records["1"].line_bytes,
            ),
        )
        self.assertEqual("chrM", records["MT"].name)

    def test_fai_reader_rejects_incomplete_and_alias_duplicate_records(self) -> None:
        fai = self.directory / "reference.fai"
        fai.write_text("1\t10\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "needs five columns"):
            read_fai(fai)
        fai.write_text("chrM\t10\t4\t5\t6\nMT\t10\t14\t5\t6\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate canonical contig 'MT'"):
            read_fai(fai)


if __name__ == "__main__":
    unittest.main()
