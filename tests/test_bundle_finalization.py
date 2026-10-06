"""Tests for final bundle eligibility and diagnostic classification."""

from __future__ import annotations

import unittest
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_finalization import classify_units, variant_class_counts  # noqa: E402

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
                "status": "zero_event_mass" if (unit, role) == ("X_UNIT", "hem_chrX") else "ok",
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
                    "allele_state": ["observed", "absent", "rejected", "observed", "rejected"],
                }
            ).write_parquet(states)
            pl.DataFrame(
                {
                    "CHROM": ["1"], "POS": [3], "REF": ["AA"], "ALT": ["A"],
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


if __name__ == "__main__":
    unittest.main()
