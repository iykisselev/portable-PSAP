"""Independent contracts for cohort-independent generation preparation."""

from __future__ import annotations

import json
import random
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_prep import (  # noqa: E402
    coverage_rows,
    effective_calibration_mask,
    read_bed_units,
    read_intervals,
    read_population_af,
    supported_scope_mask,
    subtract_intervals,
)
from scalable_generation import population_rows, role_probability_expression  # noqa: E402


class GenerationPreparationTests(unittest.TestCase):
    def test_subtraction_matches_independent_base_sets(self):
        rng = random.Random(42)
        for _ in range(100):
            regions = [tuple(sorted(rng.sample(range(40), 2))) for _ in range(5)]
            exclusions = [tuple(sorted(rng.sample(range(40), 2))) for _ in range(5)]
            expected = {p for start, end in regions for p in range(start, end)} - {p for start, end in exclusions for p in range(start, end)}
            result = subtract_intervals(regions, exclusions)
            self.assertEqual({p for start, end in result for p in range(start, end)}, expected)
            self.assertEqual(sum(end - start for start, end in result), len(expected))

    def test_subtraction_handles_many_disjoint_exclusions(self):
        exclusions = [(position, position + 1) for position in range(1, 30000, 3)]
        retained = subtract_intervals([(0, 30000)], exclusions)
        self.assertEqual(sum(end - start for start, end in retained), 20000)
        self.assertEqual(len(retained), 10001)

    def test_configured_info_tag_preserves_zero_and_rejects_other_states(self) -> None:
        observed, rejected = read_population_af(
            ROOT / "tests/null_distribution/population_af.vcf", "AF_CUSTOM"
        )
        expected = json.loads((ROOT / "tests/null_distribution/fixture_expectations.json").read_text())
        self.assertEqual(observed[("1", 10, "A", "C")], 0.125)
        self.assertEqual(observed[("1", 10, "A", "G")], 0.375)
        self.assertEqual(observed[("1", 20, "T", "C")], 0.0)
        self.assertEqual(rejected[("1", 40, "C", "T")], "non_pass_filter")
        self.assertEqual(rejected[("1", 30, "G", "A")], "missing_configured_af")
        self.assertEqual(expected["af_assignment"]["configured_tag"], "AF_CUSTOM")

    def test_union_coverage_does_not_double_count_overlapping_unit_intervals(self) -> None:
        units = read_bed_units(ROOT / "tests/null_distribution/testing_units.bed")
        coverage = read_intervals(ROOT / "tests/null_distribution/coverage.bed")
        rows = {row[0]: row for row in coverage_rows(units, coverage, 0.5)}
        self.assertEqual(rows["UNIT_OVERLAP"][2:6], (10, 7, 0.7, "eligible"))
        self.assertEqual(rows["UNIT_LOW_COVERAGE"][2:6], (10, 4, 0.4, "coverage_below_minimum"))

    def test_multiallelic_scalar_tag_is_not_split_into_alt_frequencies(self) -> None:
        _, rejected = read_population_af(
            ROOT / "tests/null_distribution/population_af.vcf", "AF_SCALAR"
        )
        self.assertEqual(rejected[("1", 50, "A", "C")], "scalar_af_on_multiallelic_record")
        self.assertEqual(rejected[("1", 50, "A", "G")], "scalar_af_on_multiallelic_record")

    def test_conflicting_duplicate_panel_allele_fails(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "duplicates.vcf"
            path.write_text(
                "##fileformat=VCFv4.3\n"
                "##INFO=<ID=CUSTOM,Number=A,Type=Float,Description=\"x\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t1\t.\tA\tC\t.\tPASS\tCUSTOM=0.1\n"
                "1\t1\t.\tA\tC\t.\tPASS\tCUSTOM=0.2\n"
            )
            with self.assertRaisesRegex(ValueError, "Duplicate/conflicting"):
                read_population_af(path, "CUSTOM")

    def test_effective_mask_intersects_coverage_with_supported_nonexcluded_scope(self) -> None:
        scope = supported_scope_mask(
            {"1": 10},
            {"1": [(2, 4), (7, 8)], "UNSUPPORTED": [(0, 5)]},
        )
        mask = effective_calibration_mask(
            {"1": [(-5, 6), (6, 15)], "2": [(0, 10)]},
            scope,
        )
        self.assertEqual(scope, {"1": [(0, 2), (4, 7), (8, 10)]})
        self.assertEqual(mask, {"1": [(0, 2), (4, 7), (8, 10)]})

    def test_x_population_rows_keep_sex_specific_frequencies_separate(self) -> None:
        rows = population_rows(
            ROOT / "tests/null_distribution/population_af.vcf",
            None,
            "X",
            0,
            10,
            {"overall": "AF_CUSTOM", "female": "AF_FEMALE", "male": "AF_MALE"},
            x_af_policy="sex_specific",
            hemizygote=True,
        )
        first = rows.filter((rows["POS"] == 1) & (rows["ALT"] == "C")).row(0, named=True)
        self.assertEqual(first["allele_state"], "observed")
        self.assertEqual(first["af_overall"], 0.2)
        self.assertEqual(first["af_female"], 0.25)
        self.assertEqual(first["af_male"], 0.4)

    def test_invalid_selected_af_is_retained_for_the_calibration_universe_to_reject(self) -> None:
        rows = population_rows(
            ROOT / "tests/null_distribution/population_af.vcf",
            None,
            "1",
            0,
            60,
            {"overall": "AF_CUSTOM"},
        )
        invalid = rows.filter((rows["POS"] == 30) & (rows["ALT"] == "A")).row(0, named=True)
        self.assertEqual(invalid["allele_state"], "invalid")
        self.assertIn("missing selected AF", invalid["state_reason"])

    def test_hem_probability_uses_haploid_af_not_a_diploid_transform(self) -> None:
        import polars as pl

        values = pl.DataFrame({"effective_af": [0.2]}).select(
            role_probability_expression("het").alias("het"),
            role_probability_expression("hom").alias("hom"),
            role_probability_expression("hem_chrX").alias("hem"),
        ).row(0, named=True)
        self.assertAlmostEqual(values["het"], 0.32)
        self.assertAlmostEqual(values["hom"], 0.04)
        self.assertAlmostEqual(values["hem"], 0.2)
        with self.assertRaisesRegex(ValueError, "Unsupported generated"):
            role_probability_expression("unknown")


if __name__ == "__main__":
    unittest.main()
