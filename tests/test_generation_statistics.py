"""Bucket-local calculation coverage."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_bindings import (  # noqa: E402
    PartitionStatisticsInputs,
    PartitionStatisticsOutputs,
    PartitionStatisticsParameters,
)
from generation_statistics import calculate_partition_statistics  # noqa: E402


class GenerationStatisticsTests(unittest.TestCase):
    def test_calculation_selects_only_the_requested_bucket_alleles(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            partition_id = "1_000000000000_000000000100"
            (directory / "partitions.json").write_text(json.dumps({
                "schema_version": 1,
                "partitions": [{
                    "partition_id": partition_id, "chromosome": "1", "start": 0,
                    "end": 100, "ordinal": 0,
                }],
            }), encoding="utf-8")
            (directory / "buckets.json").write_text(json.dumps({
                "schema_version": 2,
                "buckets": [{
                    "bucket_id": "00000", "ordinal": 0, "unit_count": 1,
                    "unit_ids": ["wanted"], "partition_ids": [partition_id],
                }],
            }), encoding="utf-8")
            membership = directory / "membership.parquet"
            scores = directory / "scores.parquet"
            states = directory / "states.parquet"
            pl.DataFrame({
                "testing_unit_id": ["wanted", "other"], "CHROM": ["1", "1"],
                "POS": [10, 20], "REF": ["A", "A"], "ALT": ["C", "G"],
                "variant_class": ["snv", "snv"],
            }).write_parquet(membership)
            pl.DataFrame({
                "CHROM": ["1", "1"], "POS": [10, 20], "REF": ["A", "A"],
                "ALT": ["C", "G"], "score": [10.0, 20.0],
                "variant_class": ["snv", "snv"],
            }).write_parquet(scores)
            pl.DataFrame({
                "CHROM": ["1", "1"], "POS": [10, 20], "REF": ["A", "A"],
                "ALT": ["C", "G"], "allele_state": ["observed", "observed"],
                "af": [0.1, 0.1], "af_overall": [0.1, 0.1],
                "af_female": [0.1, 0.1], "af_male": [0.1, 0.1],
                "state_reason": ["ok", "ok"], "af_zero_imputed_tags": [[], []],
                "variant_class": ["snv", "snv"],
            }).write_parquet(states)
            output = directory / "statistics.parquet"

            calculate_partition_statistics(
                PartitionStatisticsInputs(
                    partitions=directory / "partitions.json", buckets=directory / "buckets.json",
                    allele_states=states, allele_membership=membership, allele_scores=scores,
                    mutation_rate=None, mutation_rate_index=None,
                ),
                PartitionStatisticsOutputs(statistics=output),
                PartitionStatisticsParameters(
                    x_af_policy="overall_proxy", unobserved_policy="zero",
                    x_mutation_rate_policy="sex_averaged_proxy", score_max=70,
                    bin_width=5, roles=("het", "hom"), compression="zstd",
                ),
                partition_id,
                "00000",
            )

            statistics = pl.read_parquet(output)
            self.assertEqual(statistics["testing_unit_id"].unique().to_list(), ["wanted"])
            self.assertEqual(statistics["event_count"].sum(), 2)


if __name__ == "__main__":
    unittest.main()
