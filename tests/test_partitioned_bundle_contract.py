"""Public N-by-M partitioned-bundle component layout."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partitioned_bundle_contract import (  # noqa: E402
    PARTITIONED_NULL_BUNDLE_CONTRACT,
    lookup_components,
    lookup_filename,
)


class PartitionedBundleContractTests(unittest.TestCase):
    def test_lookup_components_are_physical_model_tables(self):
        components = lookup_components(("cds", "intronic"), hemizygote=True)
        self.assertEqual("partitioned_null_bundle_v1", PARTITIONED_NULL_BUNDLE_CONTRACT)
        self.assertEqual(10, len(components))
        self.assertEqual(
            "lookup.cds.DOM-het.tsv.gz", components[("cds", "DOM-het")]
        )
        self.assertEqual(
            "lookup.intronic.X-linked-hem.tsv.gz",
            components[("intronic", "X-linked-hem")],
        )

    def test_partition_ids_cannot_make_ambiguous_filenames(self):
        for invalid in ("a.b", "../cds", "", "has space"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    lookup_filename(invalid, "DOM-het")


if __name__ == "__main__":
    unittest.main()
