"""Partition IDs remain safe in generated public filenames."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partitioned_bundle_contract import require_partition_id


class GenerationTopologyTests(unittest.TestCase):
    def test_partition_identifier_rejects_filename_ambiguity(self):
        self.assertEqual("cds", require_partition_id("cds"))
        with self.assertRaises(ValueError):
            require_partition_id("cds.part")
