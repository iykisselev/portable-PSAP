"""Partitioned lookup filenames are explicit for every model."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partitioned_bundle_contract import lookup_components


class GenerationLookupTests(unittest.TestCase):
    def test_partition_model_lookup_layout(self):
        lookups = lookup_components(("cds", "intronic"), hemizygote=False)
        self.assertEqual(8, len(lookups))
        self.assertIn(("cds", "DOM-het"), lookups)
