"""Contracts for versioned, data-dependent generation topology manifests."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_manifests import (  # noqa: E402
    BUCKET_MANIFEST_SCHEMA_VERSION,
    PARTITION_MANIFEST_SCHEMA_VERSION,
    GenerationManifestError,
    bucket_from_manifest,
    bucket_ids,
    partition_ids,
)


class GenerationManifestTests(unittest.TestCase):
    def write(self, path: Path, payload: dict) -> Path:
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_validated_manifests_publish_canonical_dynamic_wildcards(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            partitions = self.write(
                directory / "partitions.json",
                {
                    "schema_version": PARTITION_MANIFEST_SCHEMA_VERSION,
                    "partitions": [
                        {
                            "partition_id": "1_000000000000_000000000010",
                            "chromosome": "1",
                            "start": 0,
                            "end": 10,
                            "ordinal": 0,
                        }
                    ],
                },
            )
            buckets = self.write(
                directory / "buckets.json",
                {
                    "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
                    "buckets": [
                        {
                            "bucket_id": "00000",
                            "ordinal": 0,
                            "unit_ids": ["UNIT_A", "UNIT_B"],
                        }
                    ],
                },
            )

            self.assertEqual(
                ("1_000000000000_000000000010",), partition_ids(partitions)
            )
            self.assertEqual(("00000",), bucket_ids(buckets))
            self.assertEqual(("UNIT_A", "UNIT_B"), bucket_from_manifest(buckets, "00000").unit_ids)

    def test_duplicate_or_unsorted_bucket_units_are_rejected(self) -> None:
        with TemporaryDirectory() as temporary:
            path = self.write(
                Path(temporary) / "buckets.json",
                {
                    "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
                    "buckets": [
                        {
                            "bucket_id": "00000",
                            "ordinal": 0,
                            "unit_ids": ["UNIT_B", "UNIT_B"],
                        }
                    ],
                },
            )
            with self.assertRaisesRegex(GenerationManifestError, "noncanonical"):
                bucket_ids(path)


if __name__ == "__main__":
    unittest.main()
