"""The partitioned bundle keeps SNV/indel metadata in one contract."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_contract import require_supported_contract  # noqa: E402
from partitioned_bundle_contract import (  # noqa: E402
    PARTITIONED_NULL_BUNDLE_CONTRACT,
    lookup_components,
)


class PartitionedIndelBundleContractTests(unittest.TestCase):
    def test_single_contract_supports_all_enabled_models(self) -> None:
        self.assertEqual(
            PARTITIONED_NULL_BUNDLE_CONTRACT,
            require_supported_contract(PARTITIONED_NULL_BUNDLE_CONTRACT),
        )
        components = lookup_components(("cds", "intronic"), hemizygote=True)
        self.assertEqual(10, len(components))
        self.assertEqual(
            "lookup.intronic.X-linked-hem.tsv.gz",
            components[("intronic", "X-linked-hem")],
        )


if __name__ == "__main__":
    unittest.main()
