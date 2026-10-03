"""Candidate-scoped score reads and empty partition artifact contracts."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import polars as pl
from polars.testing import assert_frame_equal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

import generation_population
from generation_partition import prepare_generation_partition
from generation_population import POPULATION_SCHEMA
from generation_scores import ALLELE_COLUMNS, SCORE_SCHEMA, score_rows


class GenerationPartitionTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.partition_id = "1_000000000000_000000000010"
        self.inputs = SimpleNamespace(**{
            name: self.directory / name for name in
            ("partitions", "units", "unit_intervals", "mask", "fasta", "fasta_fai", "score_snv", "population_vcf")
        }, score_snv_index=None, population_vcf_index=None)
        self.outputs = SimpleNamespace(**{
            name: self.directory / f"{name}.parquet"
            for name in ("allele_states", "allele_membership", "allele_scores")
        })
        self.params = SimpleNamespace(
            compression="zstd", af_tags={"overall": "AF"}, x_af_policy="overall_proxy", hemizygote=False,
        )
        self.inputs.partitions.write_text(json.dumps({
            "schema_version": 1, "partitions": [{
                "partition_id": self.partition_id, "chromosome": "1", "start": 0, "end": 10, "ordinal": 0,
            }],
        }))
        self.inputs.fasta.write_text(">1\nAAAAAAAAAA\n")
        self.inputs.fasta_fai.write_text("1\t10\t3\t10\t11\n")
        pl.DataFrame({"testing_unit_id": ["a", "b", "excluded"], "status": ["eligible", "eligible", "coverage_below_minimum"]}).write_parquet(self.inputs.units)
        pl.DataFrame({
            "testing_unit_id": ["a", "b", "excluded"], "CHROM": ["1"] * 3,
            "START": [0, 1, 8], "END": [3, 6, 10],
        }).write_parquet(self.inputs.unit_intervals)
        self.write_mask([(0, 2), (4, 10)])
        self.inputs.score_snv.write_text("".join(
            f"chr1\t{position}\tA\t{alt}\t0\t{position}.0\n"
            for position in range(1, 11) for alt in "CGT"
        ))
        self.inputs.population_vcf.write_text(
            '##fileformat=VCFv4.2\n##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            '1\t1\t.\tA\tC\t.\tPASS\tAF=0.1\n'
        )

    def write_mask(self, intervals):
        pl.DataFrame(
            [("1", start, end) for start, end in intervals],
            schema={"CHROM": pl.String, "START": pl.Int64, "END": pl.Int64}, orient="row",
        ).write_parquet(self.inputs.mask)

    def prepare(self):
        prepare_generation_partition(self.inputs, self.outputs, self.params, self.partition_id)

    def test_plain_scores_merge_clip_and_filter_intervals(self):
        scores = score_rows(self.inputs.score_snv, None, "1", 1, 6, intervals=[(4, 20), (0, 2), (1, 2)])
        self.assertEqual(scores["POS"].unique().sort().to_list(), [2, 5, 6])
        self.assertEqual(scores.height, 9)
        expected = score_rows(self.inputs.score_snv, None, "1", 1, 6).filter(pl.col("POS").is_in([2, 5, 6]))
        assert_frame_equal(scores, expected)

    def test_indexed_scores_open_once_and_fetch_only_merged_intervals(self):
        reader = Mock()
        reader.contigs = ["chr1"]
        reader.fetch.side_effect = [
            ["chr1\t2\tA\tC\t0\t2\n"] * 2,
            ["chr1\t6\tA\tC\t0\t6\n"],
        ]
        fake_pysam = Mock()
        fake_pysam.TabixFile.return_value.__enter__ = Mock(return_value=reader)
        fake_pysam.TabixFile.return_value.__exit__ = Mock(return_value=False)
        with patch.object(generation_population, "pysam", fake_pysam):
            scores = score_rows("scores.gz", "custom.tbi", "1", 1, 6, intervals=[(4, 9), (0, 2), (1, 2)])
        fake_pysam.TabixFile.assert_called_once_with("scores.gz", index="custom.tbi")
        self.assertEqual([call.args for call in reader.fetch.call_args_list], [("chr1", 1, 2), ("chr1", 4, 6)])
        # Do not conceal real duplicate source records while avoiding repeat queries.
        self.assertEqual(scores.height, 3)

    @unittest.skipUnless(generation_population.pysam, "pysam is required for real Tabix coverage")
    def test_real_tabix_and_plain_reads_agree_at_boundaries(self):
        pysam = generation_population.pysam
        compressed = str(self.directory / "scores.gz")
        pysam.tabix_compress(str(self.inputs.score_snv), compressed, force=True)
        pysam.tabix_index(compressed, seq_col=0, start_col=1, end_col=1, force=True)
        options = dict(intervals=[(4, 9), (0, 2), (1, 2)])
        expected = score_rows(self.inputs.score_snv, None, "1", 1, 6, **options)
        actual = score_rows(compressed, compressed + ".tbi", "1", 1, 6, **options)
        assert_frame_equal(actual, expected)

    def test_no_score_intervals_does_not_open_resource(self):
        for intervals in ([], [(20, 30)]):
            scores = score_rows("does-not-exist", "does-not-exist.tbi", "1", 0, 10, intervals=intervals)
            assert_frame_equal(scores, pl.DataFrame(schema=SCORE_SCHEMA))

    def test_partition_limits_scores_to_eligible_masked_units(self):
        with patch("generation_partition.score_rows", wraps=score_rows) as read_scores:
            self.prepare()
        self.assertEqual(read_scores.call_args.kwargs["intervals"], [(0, 2), (4, 6)])
        scores = pl.read_parquet(self.outputs.allele_scores)
        expected = score_rows(self.inputs.score_snv, None, "1", 0, 10).filter(pl.col("POS").is_in([1, 2, 5, 6]))
        assert_frame_equal(scores.sort(ALLELE_COLUMNS), expected.sort(ALLELE_COLUMNS))
        memberships = pl.read_parquet(self.outputs.allele_membership)
        self.assertEqual(memberships.height, 15)
        self.assertEqual(memberships.filter(pl.col("POS") == 2)["testing_unit_id"].unique().sort().to_list(), ["a", "b"])
        states = pl.read_parquet(self.outputs.allele_states)
        self.assertEqual(states.height, 12)
        self.assertEqual(states.filter(pl.col("allele_state") == "observed")["af"].to_list(), [0.1])
        self.assertEqual(states.filter(pl.col("allele_state") == "absent").height, 11)

    def test_missing_candidate_score_still_fails(self):
        lines = self.inputs.score_snv.read_text().splitlines(keepends=True)
        self.inputs.score_snv.write_text("".join(lines[1:]))
        with self.assertRaisesRegex(ValueError, "missing 1 calibration alleles"):
            self.prepare()

    def test_partition_outputs_match_full_window_score_read(self):
        self.prepare()
        optimized = {name: pl.read_parquet(path) for name, path in vars(self.outputs).items()}

        def read_full_window(*args, **kwargs):
            return score_rows(*args)

        with patch("generation_partition.score_rows", side_effect=read_full_window):
            self.prepare()
        for name, path in vars(self.outputs).items():
            with self.subTest(artifact=name):
                original = pl.read_parquet(path)
                assert_frame_equal(optimized[name].sort(original.columns), original.sort(original.columns))

    def test_nonfinite_candidate_score_still_fails(self):
        content = self.inputs.score_snv.read_text()
        self.inputs.score_snv.write_text(content.replace("\t1.0\n", "\tnan\n", 1))
        with self.assertRaisesRegex(ValueError, "non-finite score"):
            self.prepare()

    def test_duplicate_candidate_score_still_fails(self):
        content = self.inputs.score_snv.read_text()
        self.inputs.score_snv.write_text(content + content.splitlines(keepends=True)[0])
        with self.assertRaisesRegex(ValueError, "duplicate calibration alleles"):
            self.prepare()

    def test_empty_candidates_write_typed_outputs_without_resource_reads(self):
        for scenario in ("no_overlap", "all_n", "no_eligible_units"):
            with self.subTest(scenario=scenario):
                self.write_mask([(8, 10)] if scenario == "no_overlap" else [(0, 10)])
                self.inputs.fasta.write_text(">1\n" + ("N" if scenario == "all_n" else "A") * 10 + "\n")
                if scenario == "no_eligible_units":
                    pl.read_parquet(self.inputs.units).with_columns(pl.lit("coverage_below_minimum").alias("status")).write_parquet(self.inputs.units)
                with patch("generation_partition.score_rows") as scores, patch("generation_partition.population_rows") as population:
                    self.prepare()
                    scores.assert_not_called()
                    population.assert_not_called()
                assert_frame_equal(pl.read_parquet(self.outputs.allele_states), pl.DataFrame(schema=POPULATION_SCHEMA))
                assert_frame_equal(pl.read_parquet(self.outputs.allele_scores), pl.DataFrame(schema=SCORE_SCHEMA))
                assert_frame_equal(pl.read_parquet(self.outputs.allele_membership), pl.DataFrame(schema={
                    "testing_unit_id": pl.String, "CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String,
                }))


if __name__ == "__main__":
    unittest.main()
