"""Independent contracts for cohort-independent generation preparation."""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_prep import (  # noqa: E402
    retention_rows,
    effective_calibration_mask,
    read_bed_units,
    read_intervals,
    supported_scope_mask,
    subtract_intervals,
)
from generation_population import population_rows  # noqa: E402
from generation_probability import role_probability_expression  # noqa: E402
from generation_scope import prepare_generation_scope


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

    def test_scope_unions_all_exclusions_without_shrinking_retention_denominator(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            extra = directory / "extra.bed"
            extra.write_text("1\t0\t3\n1\t6\t9\n1\t20\t30\nX\t0\t4\n")
            fixture = ROOT / "tests/null_distribution"
            inputs = SimpleNamespace(
                population_vcf=fixture / "population_af.vcf", score_snv=fixture / "generation_scores.tsv",
                mutation_rate=[], exclusion=fixture / "exclusion.bed", exclusions=[extra],
                testing_units=fixture / "testing_units.bed", par_regions=fixture / "par_regions.bed",
                fasta=ROOT / "tests/smoke/data/reference.fa", fasta_fai=ROOT / "tests/smoke/data/reference.fa.fai",
            )
            outputs = SimpleNamespace(**{name: directory / name for name in (
                "partitions", "buckets", "units", "unit_intervals", "masked_unit_intervals", "mask", "resource_identities",
            )})
            params = SimpleNamespace(minimum_retained_fraction=0.5, partition_size=10, unit_bucket_size=250, compression="zstd")
            prepare_generation_scope(inputs, outputs, params)
            rows = {row["testing_unit_id"]: row for row in pl.read_parquet(outputs.units).to_dicts()}
            self.assertEqual(rows["UNIT_OVERLAP"]["supported_bases"], 10)
            self.assertEqual(rows["UNIT_OVERLAP"]["retained_bases"], 3)
            self.assertEqual(rows["UNIT_OVERLAP"]["status"], "retained_fraction_below_minimum")
            self.assertEqual(rows["UNIT_LOW_COVERAGE"]["status"], "fully_excluded")
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["supported_bases"], 8)
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["retained_fraction"], 0.5)
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["status"], "eligible")

    def test_configured_info_tag_preserves_zero_and_rejects_other_states(self) -> None:
        rows = population_rows(
            ROOT / "tests/null_distribution/population_af.vcf",
            None,
            "1",
            0,
            60,
            {"overall": "AF_CUSTOM"},
        )
        states = {
            (row["POS"], row["REF"], row["ALT"]): row
            for row in rows.to_dicts()
        }
        self.assertEqual(states[(10, "A", "C")]["af"], 0.125)
        self.assertEqual(states[(10, "A", "G")]["af"], 0.375)
        self.assertEqual(states[(20, "T", "C")]["af"], 0.0)
        self.assertEqual(states[(40, "C", "T")]["state_reason"], "non_pass_filter")
        self.assertEqual(states[(30, "G", "A")]["allele_state"], "invalid")
        self.assertIn("missing selected AF", states[(30, "G", "A")]["state_reason"])

    def test_union_retention_does_not_double_count_overlapping_unit_intervals(self) -> None:
        units = read_bed_units(ROOT / "tests/null_distribution/testing_units.bed")
        exclusions = read_intervals(ROOT / "tests/null_distribution/exclusion.bed")
        supported = supported_scope_mask({"1": 60, "X": 60}, read_intervals(ROOT / "tests/null_distribution/par_regions.bed"))
        retained = effective_calibration_mask(exclusions, supported)
        rows = {row[0]: row for row in retention_rows(units, retained, 0.5, supported)}
        self.assertEqual(rows["UNIT_OVERLAP"][2:6], (10, 7, 0.7, "eligible"))
        self.assertEqual(rows["UNIT_LOW_COVERAGE"][2:6], (10, 4, 0.4, "retained_fraction_below_minimum"))
        self.assertEqual(rows["UNIT_X_MIXED_PAR"][2:6], (8, 8, 1.0, "eligible"))

    def test_multiallelic_scalar_tag_is_not_split_into_alt_frequencies(self) -> None:
        rows = population_rows(
            ROOT / "tests/null_distribution/population_af.vcf",
            None,
            "1",
            0,
            60,
            {"overall": "AF_SCALAR"},
        )
        invalid = rows.filter((pl.col("POS") == 50) & pl.col("ALT").is_in(["C", "G"]))
        self.assertEqual(invalid["allele_state"].unique().to_list(), ["invalid"])
        self.assertTrue(
            invalid["state_reason"].str.contains("scalar selected AF on a multiallelic record").all()
        )

    def test_effective_mask_subtracts_exclusions_from_supported_nonpar_scope(self) -> None:
        scope = supported_scope_mask(
            {"1": 10},
            {"1": [(2, 4), (7, 8)], "UNSUPPORTED": [(0, 5)]},
        )
        mask = effective_calibration_mask(
            {"1": [(1, 5), (4, 6)], "2": [(0, 10)]},
            scope,
        )
        self.assertEqual(scope, {"1": [(0, 2), (4, 7), (8, 10)]})
        self.assertEqual(mask, {"1": [(0, 1), (6, 7), (8, 10)]})

    def test_retention_threshold_counts_exclusions_only_in_numerator(self) -> None:
        units = {"half": [("X", 0, 12)], "low": [("1", 0, 10)],
                 "removed": [("1", 20, 30)], "par": [("X", 4, 6)]}
        supported = supported_scope_mask({"1": 40, "X": 12}, {"X": [(4, 6)]})
        retained = effective_calibration_mask({"X": [(0, 4), (6, 7)], "1": [(4, 10), (20, 30)]}, supported)
        rows = {row[0]: row for row in retention_rows(units, retained, 0.5, supported)}
        self.assertEqual(rows["half"][2:], (10, 5, 0.5, "eligible"))
        self.assertEqual(rows["low"][2:], (10, 4, 0.4, "retained_fraction_below_minimum"))
        self.assertEqual(rows["removed"][2:], (10, 0, 0.0, "fully_excluded"))
        self.assertEqual(rows["par"][2:], (0, 0, 0.0, "no_calibration_scope"))

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
