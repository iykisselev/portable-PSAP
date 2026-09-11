"""Tests for final bundle eligibility and diagnostic classification."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workflow/scripts/generate/finalize_bundle.py"
SPEC = importlib.util.spec_from_file_location("finalize_bundle", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

class BundleFinalizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.units = [
            {"testing_unit_id": "ELIGIBLE", "status": "eligible"},
            {"testing_unit_id": "LOW", "status": "coverage_below_minimum"},
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
        eligible, excluded = MODULE.classify_units(
            self.units, self.diagnostics("zero_event_mass", "zero_event_mass")
        )
        self.assertEqual(eligible, [])
        self.assertEqual(
            excluded,
            [
                ("ELIGIBLE", "zero_event_mass_het_hom"),
                ("LOW", "coverage_below_minimum"),
            ],
        )

    def test_numerical_failure_aborts_finalization(self) -> None:
        with self.assertRaisesRegex(ValueError, "Numerical distribution failure"):
            MODULE.classify_units(self.units, self.diagnostics("numerical_failure"))

    def test_missing_diagnostic_aborts_finalization(self) -> None:
        with self.assertRaisesRegex(ValueError, "Missing distribution diagnostic"):
            MODULE.classify_units(self.units, self.diagnostics()[:1])

    def test_duplicate_diagnostic_aborts_finalization(self) -> None:
        diagnostics = self.diagnostics()
        with self.assertRaisesRegex(ValueError, "Duplicate distribution diagnostic"):
            MODULE.classify_units(self.units, diagnostics + diagnostics[:1])

    def test_malformed_event_mass_aborts_finalization(self) -> None:
        diagnostics = self.diagnostics()
        diagnostics[0]["event_mass"] = "not-a-number"
        with self.assertRaisesRegex(ValueError, "Invalid event mass"):
            MODULE.classify_units(self.units, diagnostics)

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
        eligible, excluded = MODULE.classify_units(units, diagnostics, hemizygote=True)
        self.assertEqual(eligible, ["AUTO"])
        self.assertEqual(excluded, [("X_UNIT", "zero_event_mass_hem_chrX")])


if __name__ == "__main__":
    unittest.main()
