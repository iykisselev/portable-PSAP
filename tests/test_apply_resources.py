"""Apply resources are resolved from the certified partitioned bundle."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from partitioned_bundle_contract import PARTITIONED_NULL_BUNDLE_CONTRACT


class ApplyResourceTests(unittest.TestCase):
    def test_only_partitioned_bundle_contract_is_supported(self):
        self.assertEqual("partitioned_null_bundle_v1", PARTITIONED_NULL_BUNDLE_CONTRACT)
