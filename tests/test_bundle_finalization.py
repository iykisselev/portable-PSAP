"""Tests for final bundle eligibility and diagnostic classification."""

from __future__ import annotations

import csv
import gzip
import unittest
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_finalization import (  # noqa: E402
    classify_units,
    variant_class_counts,
    write_finalization_allele_artifacts,
)


def allele_outputs(directory: Path) -> SimpleNamespace:
    return SimpleNamespace(
        excluded_alleles=str(directory / "excluded_alleles.tsv.gz"),
        candidate_audit=str(directory / "candidate_audit.tsv.gz"),
        population_af_zero_imputations=str(directory / "zero_imputations.tsv.gz"),
    )


def gzip_rows(path: str) -> list[list[str]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter="\t"))[1:]


class BundleFinalizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.units = [
            {"testing_unit_id": "ELIGIBLE", "status": "eligible"},
            {"testing_unit_id": "LOW", "status": "retained_fraction_below_minimum"},
        ]

    def diagnostics(self, het="ok", hom="ok"):
        return [
            {
                "testing_unit_id": "ELIGIBLE",
                "role": "het",
                "status": het,
                "event_mass": "0" if het == "zero_event_mass" else "0.25",
            },
            {
                "testing_unit_id": "ELIGIBLE",
                "role": "hom",
                "status": hom,
                "event_mass": "0" if hom == "zero_event_mass" else "0.25",
            },
        ]

    def test_only_true_zero_mass_is_reported_as_biological_exclusion(self) -> None:
        eligible, excluded = classify_units(
            self.units, self.diagnostics("zero_event_mass", "zero_event_mass")
        )
        self.assertEqual(eligible, [])
        self.assertEqual(
            excluded,
            [
                ("ELIGIBLE", "zero_event_mass_het_hom"),
                ("LOW", "retained_fraction_below_minimum"),
            ],
        )

    def test_numerical_failure_aborts_finalization(self) -> None:
        with self.assertRaisesRegex(ValueError, "Numerical distribution failure"):
            classify_units(self.units, self.diagnostics("numerical_failure"))

    def test_missing_diagnostic_aborts_finalization(self) -> None:
        with self.assertRaisesRegex(ValueError, "Missing distribution diagnostic"):
            classify_units(self.units, self.diagnostics()[:1])

    def test_duplicate_diagnostic_aborts_finalization(self) -> None:
        diagnostics = self.diagnostics()
        with self.assertRaisesRegex(ValueError, "Duplicate distribution diagnostic"):
            classify_units(self.units, diagnostics + diagnostics[:1])

    def test_malformed_event_mass_aborts_finalization(self) -> None:
        diagnostics = self.diagnostics()
        diagnostics[0]["event_mass"] = "not-a-number"
        with self.assertRaisesRegex(ValueError, "Invalid event mass"):
            classify_units(self.units, diagnostics)

    def test_hem_role_is_required_only_for_chromosome_x_units(self) -> None:
        units = [
            {"testing_unit_id": "AUTO", "CHROM": "1", "status": "eligible"},
            {"testing_unit_id": "X_UNIT", "CHROM": "X", "status": "eligible"},
        ]
        diagnostics = [
            {
                "testing_unit_id": unit,
                "role": role,
                "status": "zero_event_mass"
                if (unit, role) == ("X_UNIT", "hem_chrX")
                else "ok",
                "event_mass": "0" if (unit, role) == ("X_UNIT", "hem_chrX") else "0.25",
            }
            for unit in ("AUTO", "X_UNIT")
            for role in ("het", "hom", "hem_chrX")
        ]
        eligible, excluded = classify_units(units, diagnostics, hemizygote=True)
        self.assertEqual(eligible, ["AUTO"])
        self.assertEqual(excluded, [("X_UNIT", "zero_event_mass_hem_chrX")])

    def test_variant_class_counts_include_structural_indel_rejections(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            states = directory / "states.parquet"
            audit = directory / "audit.parquet"
            pl.DataFrame(
                {
                    "variant_class": ["snv", "snv", "snv", "indel", "indel"],
                    "allele_state": [
                        "observed",
                        "absent",
                        "rejected",
                        "observed",
                        "rejected",
                    ],
                }
            ).write_parquet(states)
            pl.DataFrame(
                {
                    "CHROM": ["1"],
                    "POS": [3],
                    "REF": ["AA"],
                    "ALT": ["A"],
                    "rejection_reason": ["crosses_mask_boundary"],
                }
            ).write_parquet(audit)

            self.assertEqual(
                {
                    "snv": {"retained": 2, "rejected": 1, "observed": 1, "absent": 1},
                    "indel": {"retained": 1, "rejected": 2, "observed": 1, "absent": 0},
                },
                variant_class_counts([states], [audit]),
            )

    def test_finalization_rejects_unsupported_population_state(self) -> None:
        with TemporaryDirectory() as temporary:
            states = Path(temporary) / "states.parquet"
            audit = Path(temporary) / "audit.parquet"
            pl.DataFrame(
                {
                    "CHROM": ["1"],
                    "POS": [1],
                    "REF": ["A"],
                    "ALT": ["C"],
                    "allele_state": ["invalid"],
                    "state_reason": ["bad_af"],
                    "af_zero_imputed_tags": [None],
                    "variant_class": ["snv"],
                }
            ).write_parquet(states)
            pl.DataFrame(
                schema={
                    "CHROM": pl.String,
                    "POS": pl.Int64,
                    "REF": pl.String,
                    "ALT": pl.String,
                    "rejection_reason": pl.String,
                }
            ).write_parquet(audit)
            with self.assertRaisesRegex(ValueError, "Unsupported allele state"):
                write_finalization_allele_artifacts(
                    [states],
                    [audit],
                    allele_outputs(Path(temporary)),
                    include_variant_class_counts=True,
                )

    def test_finalization_preserves_partition_order_across_sequential_state_scans(
        self,
    ) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            first, second = directory / "first.parquet", directory / "second.parquet"
            audit = directory / "audit.parquet"
            schema = {
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
                "allele_state": pl.String,
                "state_reason": pl.String,
                "af_zero_imputed_tags": pl.String,
                "variant_class": pl.String,
            }
            pl.DataFrame(
                [("1", 2, "A", "C", "rejected", "first", None, "snv")],
                schema=schema,
                orient="row",
            ).write_parquet(first)
            pl.DataFrame(
                [("1", 5, "G", "T", "rejected", "second", "AF", "snv")],
                schema=schema,
                orient="row",
            ).write_parquet(second)
            pl.DataFrame(
                schema={
                    "CHROM": pl.String,
                    "POS": pl.Int64,
                    "REF": pl.String,
                    "ALT": pl.String,
                    "rejection_reason": pl.String,
                }
            ).write_parquet(audit)

            outputs = allele_outputs(directory)
            write_finalization_allele_artifacts(
                [first, second], [audit], outputs, include_variant_class_counts=False
            )

            self.assertEqual(
                [["1", "2", "A", "C", "first"], ["1", "5", "G", "T", "second"]],
                gzip_rows(outputs.excluded_alleles),
            )
            self.assertEqual(
                [["1", "5", "G", "T", "AF"]],
                gzip_rows(outputs.population_af_zero_imputations),
            )

    def test_finalization_collects_public_allele_rows_and_v3_counts_with_native_state_scans(
        self,
    ) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            states = directory / "states.parquet"
            audit = directory / "audit.parquet"
            pl.DataFrame(
                {
                    "CHROM": ["1", "1", "X"],
                    "POS": [1, 2, 3],
                    "REF": ["A", "C", "G"],
                    "ALT": ["T", "G", "GA"],
                    "allele_state": ["observed", "rejected", "absent"],
                    "state_reason": ["present", "filtered", "not_present"],
                    "af_zero_imputed_tags": [None, "AF", "AF_FEMALE"],
                    "variant_class": ["snv", "snv", "indel"],
                }
            ).write_parquet(states)
            pl.DataFrame(
                {
                    "CHROM": ["X"],
                    "POS": [4],
                    "REF": ["AA"],
                    "ALT": ["A"],
                    "rejection_reason": ["crosses_mask_boundary"],
                }
            ).write_parquet(audit)

            outputs = allele_outputs(directory)
            _, _, _, _, _, counts = write_finalization_allele_artifacts(
                [states], [audit], outputs, include_variant_class_counts=True
            )

            self.assertEqual(
                [
                    ["1", "2", "C", "G", "filtered"],
                    ["X", "4", "AA", "A", "crosses_mask_boundary"],
                ],
                gzip_rows(outputs.excluded_alleles),
            )
            self.assertEqual(
                [["1", "2", "C", "G", "AF"], ["X", "3", "G", "GA", "AF_FEMALE"]],
                gzip_rows(outputs.population_af_zero_imputations),
            )
            self.assertEqual(
                [["X", "4", "AA", "A", "crosses_mask_boundary"]],
                gzip_rows(outputs.candidate_audit),
            )
            self.assertEqual(
                {
                    "snv": {"retained": 1, "rejected": 1, "observed": 1, "absent": 0},
                    "indel": {"retained": 1, "rejected": 1, "observed": 0, "absent": 1},
                },
                counts,
            )


if __name__ == "__main__":
    unittest.main()
