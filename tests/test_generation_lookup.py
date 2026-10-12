"""Bucket diagnostics retain independently auditable SNV and indel event counts."""

from __future__ import annotations

import csv
import json
from math import log1p
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_bindings import (  # noqa: E402
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)
from generation_lookup import reduce_partition_statistics  # noqa: E402


class GenerationLookupTests(unittest.TestCase):
    def test_bucket_diagnostics_split_nonzero_events_by_variant_class(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            buckets = directory / "buckets.json"
            statistics = directory / "statistics.parquet"
            diagnostics = directory / "diagnostics.tsv"
            buckets.write_text(
                json.dumps(
                    {
                        "schema_version": 2,
                        "buckets": [
                            {
                                "bucket_id": "00000",
                                "ordinal": 0,
                                "unit_count": 1,
                                "unit_ids": ["UNIT"],
                                "partition_ids": ["1_000000000000_000000000010"],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            pl.DataFrame(
                [
                    ("UNIT", "het", 1, 0, log1p(-0.1), 1, 1, 0),
                    ("UNIT", "het", 2, 0, log1p(-0.2), 1, 0, 1),
                    ("UNIT", "hom", 1, 0, log1p(-0.01), 1, 1, 0),
                    ("UNIT", "hom", 2, 0, log1p(-0.04), 1, 0, 1),
                ],
                schema=[
                    "testing_unit_id",
                    "role",
                    "bin",
                    "certain_count",
                    "log_none_sum",
                    "event_count",
                    "snv_event_count",
                    "indel_event_count",
                ],
                orient="row",
            ).write_parquet(statistics)

            reduce_partition_statistics(
                LookupReductionInputs(
                    buckets=buckets, partition_statistics=(statistics,)
                ),
                LookupReductionOutputs(
                    lookup_het=directory / "het.tsv",
                    lookup_hom=directory / "hom.tsv",
                    lookup_hem_chrX=None,
                    diagnostics=diagnostics,
                ),
                LookupReductionParameters(bin_count=2, roles=("het", "hom")),
                "00000",
            )

            with diagnostics.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [
                    "testing_unit_id",
                    "role",
                    "event_mass",
                    "status",
                    "event_count",
                    "snv_event_count",
                    "indel_event_count",
                ],
                list(rows[0]),
            )
            self.assertEqual(
                [
                    (
                        row["role"],
                        row["event_count"],
                        row["snv_event_count"],
                        row["indel_event_count"],
                    )
                    for row in rows
                ],
                [("het", "2", "1", "1"), ("hom", "2", "1", "1")],
            )


if __name__ == "__main__":
    unittest.main()
