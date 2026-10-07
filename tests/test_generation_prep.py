"""Independent contracts for cohort-independent generation preparation."""

from __future__ import annotations

import json
import random
import sys
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_prep import (  # noqa: E402
    overlap_sorted,
    retention_rows,
    effective_calibration_mask,
    read_bed_units,
    read_intervals,
    supported_scope_mask,
    subtract_intervals,
)
from generation_population import population_rows, validate_population_af_headers  # noqa: E402
from generation_probability import role_probability_expression, unit_bucket_frame  # noqa: E402
from generation_bindings import (  # noqa: E402
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)
from generation_scope import prepare_generation_scope


class GenerationPreparationTests(unittest.TestCase):
    def test_unit_buckets_follow_genomic_order_not_unit_id(self):
        units = pl.DataFrame({
            "testing_unit_id": ["ZETA", "ALPHA", "MID", "excluded"],
            "status": ["eligible", "eligible", "eligible", "fully_excluded"],
        })
        intervals = pl.DataFrame({
            "testing_unit_id": ["ZETA", "ALPHA", "MID", "excluded"],
            "CHROM": ["1", "2", "1", "1"],
            "START": [50, 1, 10, 0],
            "END": [60, 2, 20, 1],
        })
        self.assertEqual(
            unit_bucket_frame(units, 2, intervals).rows(),
            [("MID", 0), ("ZETA", 0), ("ALPHA", 1)],
        )

    def test_sorted_overlap_finds_only_overlapping_mask_regions(self):
        intervals = [(35, 48), (5, 17), (14, 23)]
        mask = [(0, 10), (15, 20), (30, 40), (45, 50)]
        self.assertEqual(
            overlap_sorted(intervals, mask, [start for start, _ in mask]),
            [(5, 10), (15, 20), (35, 40), (45, 48)],
        )

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
            inputs = GenerationScopeInputs(
                population_vcf=fixture / "population_af.vcf",
                population_vcf_index=None,
                score_snv=fixture / "generation_scores.tsv",
                score_snv_index=None,
                mutation_rate=None,
                mutation_rate_index=None,
                exclusion=fixture / "exclusion.bed",
                exclusions={"extra": extra},
                testing_units=fixture / "testing_units.bed",
                par_regions=fixture / "par_regions.bed",
                fasta=ROOT / "tests/smoke/data/reference.fa",
                fasta_fai=ROOT / "tests/smoke/data/reference.fa.fai",
            )
            outputs = GenerationScopeOutputs(
                partitions=directory / "partitions",
                buckets=directory / "buckets",
                units=directory / "units",
                unit_intervals=directory / "unit_intervals",
                masked_unit_intervals=directory / "masked_unit_intervals",
                mask=directory / "mask",
                resource_identities=directory / "resource_identities",
            )
            params = GenerationScopeParameters(
                0.5, 10, 250, "zstd", {"overall": "AF_CUSTOM"}
            )
            prepare_generation_scope(inputs, outputs, params)
            rows = {row["testing_unit_id"]: row for row in pl.read_parquet(outputs.units).to_dicts()}
            self.assertEqual(rows["UNIT_OVERLAP"]["supported_bases"], 10)
            self.assertEqual(rows["UNIT_OVERLAP"]["retained_bases"], 3)
            self.assertEqual(rows["UNIT_OVERLAP"]["status"], "retained_fraction_below_minimum")
            self.assertEqual(rows["UNIT_LOW_COVERAGE"]["status"], "fully_excluded")
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["supported_bases"], 8)
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["retained_fraction"], 0.5)
            self.assertEqual(rows["UNIT_X_MIXED_PAR"]["status"], "eligible")
            partitions = json.loads(outputs.partitions.read_text(encoding="utf-8"))["partitions"]
            self.assertEqual(
                ["X_000000000000_000000000010"],
                [partition["partition_id"] for partition in partitions],
            )

    def test_indel_scope_wiring_adds_only_the_required_owner_halo_partition(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fasta = directory / "reference.fa"
            fasta_index = directory / "reference.fa.fai"
            units = directory / "units.bed"
            population = directory / "population.vcf"
            scores = directory / "scores.tsv"
            indels = directory / "indels.tsv"
            par_regions = directory / "par.bed"
            fasta.write_text(">1\n" + "A" * 30 + "\n", encoding="utf-8")
            fasta_index.write_text("1\t30\t0\t0\t0\n", encoding="utf-8")
            units.write_text("1\t10\t12\tUNIT\n", encoding="utf-8")
            population.write_text(
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                encoding="utf-8",
            )
            scores.write_text("1\t11\tA\tC\t0\t20\n", encoding="utf-8")
            indels.write_text("1\t9\tA\tAT\t0\t20\n", encoding="utf-8")
            par_regions.write_text("", encoding="utf-8")
            inputs = GenerationScopeInputs(
                population_vcf=population, population_vcf_index=None,
                score_snv=scores, score_snv_index=None,
                exclusion=None, exclusions={}, testing_units=units, par_regions=par_regions,
                fasta=fasta, fasta_fai=fasta_index, mutation_rate=None, mutation_rate_index=None,
                score_indel=indels, score_indel_index=None,
            )
            outputs = GenerationScopeOutputs(
                partitions=directory / "partitions.json", buckets=directory / "buckets.json",
                units=directory / "units.parquet", unit_intervals=directory / "unit_intervals.parquet",
                masked_unit_intervals=directory / "masked.parquet", mask=directory / "mask.parquet",
                resource_identities=directory / "resources.json",
            )
            prepare_generation_scope(
                inputs, outputs,
                GenerationScopeParameters(0.5, 10, 10, "zstd", {"overall": "AF"}, True, 3),
            )
            partitions = json.loads(outputs.partitions.read_text(encoding="utf-8"))["partitions"]
            buckets = json.loads(outputs.buckets.read_text(encoding="utf-8"))["buckets"]
        self.assertEqual(
            ["1_000000000000_000000000010", "1_000000000010_000000000020"],
            [entry["partition_id"] for entry in partitions],
        )
        self.assertEqual(
            ["1_000000000000_000000000010", "1_000000000010_000000000020"],
            buckets[0]["partition_ids"],
        )

    def test_configured_info_tag_preserves_zero_and_rejects_other_states(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
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
        self.assertEqual(states[(50, "A", "C")]["af"], 0.0)
        self.assertEqual(states[(50, "A", "G")]["af_zero_imputed_tags"], "AF_CUSTOM")
        self.assertEqual(states[(50, "A", "C")]["state_reason"], "selected_af_zero_imputed")
        self.assertTrue(any("using AF=0" in str(warning.message) for warning in caught))

    def test_configured_af_header_is_required_before_partition_work(self) -> None:
        with self.assertRaisesRegex(ValueError, "no INFO header for configured tag 'AF_MISSING'"):
            validate_population_af_headers(
                ROOT / "tests/null_distribution/population_af.vcf",
                {"overall": "AF_CUSTOM", "female": "AF_MISSING"},
            )

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
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

    def test_multiallelic_selected_af_sum_must_not_exceed_one(self) -> None:
        with TemporaryDirectory() as temporary:
            population = Path(temporary) / "population.vcf"
            population.write_text(
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t2\t.\tA\tC,G\t.\tPASS\tAF=0.75,0.75\n",
                encoding="utf-8",
            )
            rows = population_rows(population, None, "1", 0, 10, {"overall": "AF"})
        self.assertEqual(rows["allele_state"].unique().to_list(), ["invalid"])
        self.assertTrue(rows["state_reason"].str.contains("AF sum exceeds 1").all())

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
