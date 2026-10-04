"""Compatibility contracts for PED and sample-sheet metadata readers."""

from __future__ import annotations

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from validation import read_ped, read_samples_tsv  # noqa: E402


class SampleMetadataCharacterizationTests(unittest.TestCase):
    def read(self, reader, path: Path) -> tuple[dict, list[str], list[str]]:
        errors: list[str] = []
        warnings: list[str] = []
        return reader(path, errors, warnings), errors, warnings

    def test_common_aggregation_keeps_format_specific_diagnostics(self) -> None:
        ped_rows = (
            "\tbad/path\t0\t0\t3\t0\n"
            "F1\taffected\t0\t0\t3\t2\n"
            "F1\tunaffected\t0\t0\t2\t1\n"
            "F1\tunaffected\t0\t0\t2\t1\n"
            "F1\t\t0\t0\t1\t2\n"
        )
        sample_rows = (
            "\tbad/path\t3\t0\n"
            "F1\taffected\t3\t2\n"
            "F1\tunaffected\t2\t1\n"
            "F1\tunaffected\t2\t1\n"
            "F1\t\t1\t2\n"
        )
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            ped = directory / "cohort.ped"
            samples = directory / "samples.tsv"
            ped.write_text(ped_rows, encoding="utf-8")
            samples.write_text(
                "family_id\tpid\tsex\tphenotype\n" + sample_rows,
                encoding="utf-8",
            )
            ped_result, ped_errors, ped_warnings = self.read(read_ped, ped)
            sample_result, sample_errors, sample_warnings = self.read(
                read_samples_tsv, samples
            )

        self.assertEqual(ped_result, sample_result)
        self.assertEqual(
            ped_errors,
            [
                "PED line 1 has an empty family ID",
                "PED line 1 sample ID contains a path separator: 'bad/path'",
                "PED line 5 has an empty sample ID",
                "Duplicate sample IDs found in PED column 2: unaffected",
            ],
        )
        self.assertEqual(
            sample_errors,
            [
                "Sample sheet line 2 has an empty family ID",
                "Sample sheet line 2 sample ID contains a path separator: 'bad/path'",
                "Sample sheet line 6 has an empty sample ID",
                "Duplicate sample IDs found in sample sheet pid column: unaffected",
            ],
        )
        self.assertEqual(
            ped_warnings,
            [
                "Excluded 1 PED sample(s) with missing or unsupported phenotype values "
                "from all processing; observed values: 0; sample IDs: bad/path",
                "PED sex values other than 1 and 2 will not receive sex-specific "
                "chromosome-X recoding: 3",
            ],
        )
        self.assertEqual(
            sample_warnings,
            [
                "Excluded 1 sample-sheet sample(s) with missing or unsupported "
                "phenotype values from all processing; observed values: 0; sample IDs: "
                "bad/path",
                "Sample-sheet sex values other than 1 and 2 will not receive sex-specific "
                "chromosome-X recoding: 3",
            ],
        )

    def test_excluded_sample_warning_keeps_ten_identifier_examples(self) -> None:
        ped_rows = "F\tkept\t0\t0\t1\t1\n" + "".join(
            f"F\tunknown{index}\t0\t0\t2\t0\n" for index in range(11)
        )
        sample_rows = "F\tkept\t1\t1\n" + "".join(
            f"F\tunknown{index}\t2\t0\n" for index in range(11)
        )
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            ped = directory / "cohort.ped"
            samples = directory / "samples.tsv"
            ped.write_text(ped_rows, encoding="utf-8")
            samples.write_text(
                "family_id\tpid\tsex\tphenotype\n" + sample_rows,
                encoding="utf-8",
            )
            _, _, ped_warnings = self.read(read_ped, ped)
            _, _, sample_warnings = self.read(read_samples_tsv, samples)

        examples = ", ".join(f"unknown{index}" for index in range(10))
        self.assertEqual(
            ped_warnings,
            [
                "Excluded 11 PED sample(s) with missing or unsupported phenotype values "
                f"from all processing; observed values: 0; sample IDs: {examples} (and 1 more)"
            ],
        )
        self.assertEqual(
            sample_warnings,
            [
                "Excluded 11 sample-sheet sample(s) with missing or unsupported phenotype "
                f"values from all processing; observed values: 0; sample IDs: {examples} (and 1 more)"
            ],
        )

    def test_format_specific_malformed_rows_and_line_numbers_are_preserved(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            short_ped = directory / "short.ped"
            short_ped.write_text("F\tsample\n", encoding="utf-8")
            _, short_errors, short_warnings = self.read(read_ped, short_ped)
            self.assertEqual(
                short_errors,
                ["PED line 1 has 2 columns; at least 6 are required"],
            )
            self.assertEqual(short_warnings, [])

            missing_header = directory / "missing.tsv"
            missing_header.write_text("family_id\tpid\tsex\n", encoding="utf-8")
            _, header_errors, header_warnings = self.read(
                read_samples_tsv, missing_header
            )
            self.assertEqual(
                header_errors,
                ["Sample sheet is missing required column(s): phenotype"],
            )
            self.assertEqual(header_warnings, [])

            ped = directory / "line-numbers.ped"
            ped.write_text("\n\tsample\t0\t0\t1\t2\n", encoding="utf-8")
            samples = directory / "line-numbers.tsv"
            samples.write_text(
                "family_id\tpid\tsex\tphenotype\n\n\tsample\t1\t2\n",
                encoding="utf-8",
            )
            _, ped_errors, _ = self.read(read_ped, ped)
            _, sample_errors, _ = self.read(read_samples_tsv, samples)
            self.assertEqual(ped_errors, ["PED line 2 has an empty family ID"])
            self.assertEqual(
                sample_errors, ["Sample sheet line 2 has an empty family ID"]
            )

    def test_input_decoding_remains_format_specific(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            ped = directory / "encoding.ped"
            ped.write_bytes(b"F\xff\tsample\t0\t0\t1\t2\n")
            result, errors, warnings = self.read(read_ped, ped)
            self.assertEqual(errors, [])
            self.assertEqual(
                warnings,
                [
                    "PED file contains no unaffected individuals (phenotype code 1); "
                    "the cohort report will use validation=no_controls"
                ],
            )
            self.assertEqual(result["manifest"][0]["family_id"], "F\ufffd")

            samples = directory / "encoding.tsv"
            samples.write_bytes(
                b"family_id\tpid\tsex\tphenotype\nF\xff\tsample\t1\t2\n"
            )
            with self.assertRaises(UnicodeDecodeError):
                self.read(read_samples_tsv, samples)


if __name__ == "__main__":
    unittest.main()
