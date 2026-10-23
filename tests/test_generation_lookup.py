"""Partitioned lookup filenames are explicit for every model."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

import polars as pl

from partitioned_bundle_contract import lookup_components, lookup_filename
from partitioned_lookup import partition_model_lookup_rows


class GenerationLookupTests(unittest.TestCase):
    def test_partition_model_lookup_layout(self):
        lookups = lookup_components(("cds", "intronic"), hemizygote=False)
        self.assertEqual(8, len(lookups))
        self.assertIn(("cds", "DOM-het"), lookups)

    def test_partition_model_tables_are_physical_and_unambiguous(self):
        lookups = lookup_components(("cds", "intronic"), hemizygote=True)
        self.assertEqual(
            "lookup.intronic.X-linked-hem.tsv.gz", lookups[("intronic", "X-linked-hem")]
        )
        for invalid in ("a.b", "../cds", "", "has space"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    lookup_filename(invalid, "DOM-het")

    def test_lookup_reduction_emits_all_models_per_partition(self):
        statistics = pl.DataFrame(
            {
                "testing_unit_id": ["AUTO", "XGENE"],
                "partition_id": ["cds", "cds"],
                "role": ["het", "het"],
                "bin": [1, 1],
                "certain_count": [0, 0],
                "log_none_sum": [-0.1, -0.2],
                "event_count": [1, 1],
                "snv_event_count": [1, 1],
                "indel_event_count": [0, 0],
            }
        )
        tables, _ = partition_model_lookup_rows(
            statistics,
            [
                {"testing_unit_id": "AUTO", "partition_id": "cds", "CHROM": "1"},
                {"testing_unit_id": "XGENE", "partition_id": "cds", "CHROM": "X"},
            ],
            hemizygote=True,
            bin_count=2,
        )
        self.assertEqual(5, len(tables))
        self.assertEqual(["AUTO"], [row[0] for row in tables[("cds", "DOM-het")]])
        self.assertEqual(["XGENE"], [row[0] for row in tables[("cds", "X-linked-het")]])
