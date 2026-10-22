"""Typed generation bindings expose partition/model outputs only."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_bindings import LookupReductionOutputs


class GenerationAdapterBindingTests(unittest.TestCase):
    def test_lookup_reduction_has_only_partition_model_outputs(self):
        output = LookupReductionOutputs(("cds.DOM.tsv",), "diagnostics.tsv")
        self.assertEqual(("cds.DOM.tsv",), output.partition_model_lookups)
