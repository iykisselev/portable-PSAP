"""Focused checks for atomic generation artifact publication helpers."""

from __future__ import annotations

import gzip
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_contract import BundleContractError, CURRENT_GENERATED_BUNDLE_CONTRACT, upgrade_metadata  # noqa: E402
from io_contracts import component_inventory, write_gzip_tsv_atomic  # noqa: E402


class IOContractTests(unittest.TestCase):
    def test_gzip_component_is_deterministic_and_has_a_complete_inventory(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "component.tsv.gz"
            count, first = write_gzip_tsv_atomic(path, ["id", "value"], [("a", 1), ("b", 2)])
            first_bytes = path.read_bytes()
            _, second = write_gzip_tsv_atomic(path, ["id", "value"], [("a", 1), ("b", 2)])
            self.assertEqual(count, 2)
            self.assertEqual(first_bytes, path.read_bytes())
            self.assertEqual(first, second)
            self.assertEqual(first, component_inventory(path))
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "id\tvalue\na\t1\nb\t2\n")

    def test_contract_upgrade_is_explicitly_not_a_silent_v1_conversion(self) -> None:
        self.assertEqual(
            upgrade_metadata({"contract": CURRENT_GENERATED_BUNDLE_CONTRACT}),
            {"contract": CURRENT_GENERATED_BUNDLE_CONTRACT},
        )
        with self.assertRaises(BundleContractError):
            upgrade_metadata({"contract": "generated_bundle_v1"})


if __name__ == "__main__":
    unittest.main()
