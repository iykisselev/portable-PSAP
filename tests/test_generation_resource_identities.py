"""Deferred generation-resource provenance coverage."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_bindings import GenerationResourceIdentityOutputs  # noqa: E402
from generation_resource_identities import hash_generation_resources  # noqa: E402
from io_contracts import sha256_file  # noqa: E402


class GenerationResourceIdentityTests(unittest.TestCase):
    def test_hashing_publishes_the_same_resource_inventory_outside_the_checkpoint(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "resource.txt"
            output = directory / "resource_identities.json"
            source.write_text("generation resource\n", encoding="utf-8")

            hash_generation_resources(
                {"score_snv": source},
                GenerationResourceIdentityOutputs(resource_identities=output),
            )

            record = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(str(source), record["resources"]["score_snv"]["path"])
            self.assertEqual(sha256_file(source), record["resources"]["score_snv"]["sha256"])


if __name__ == "__main__":
    unittest.main()
