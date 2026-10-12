"""Contracts for versioned, data-dependent generation topology manifests."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

import generation_manifests  # noqa: E402
from generation_manifests import (  # noqa: E402
    BUCKET_MANIFEST_SCHEMA_VERSION,
    PARTITION_MANIFEST_SCHEMA_VERSION,
    GenerationManifestError,
    PartitionManifestEntry,
    bucket_from_manifest,
    bucket_ids,
    partition_ids,
    read_partition_manifest,
)
from generation_partitions import (  # noqa: E402
    GenomePartition,
    build_partitions,
    indel_candidate_owner_intervals,
    read_genome_partitions,
)
from generation_scope import sparse_bucket_partitions  # noqa: E402


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
                            "partition_ids": ["1_000000000000_000000000010"],
                        }
                    ],
                },
            )

            self.assertEqual(
                ("1_000000000000_000000000010",), partition_ids(partitions)
            )
            self.assertIsInstance(
                read_partition_manifest(partitions)[0], PartitionManifestEntry
            )
            self.assertIsInstance(
                read_genome_partitions(partitions)[0], GenomePartition
            )
            self.assertEqual(("00000",), bucket_ids(buckets))
            self.assertEqual(
                ("UNIT_A", "UNIT_B"), bucket_from_manifest(buckets, "00000").unit_ids
            )
            self.assertEqual(
                ("1_000000000000_000000000010",),
                bucket_from_manifest(buckets, "00000").partition_ids,
            )

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
                            "partition_ids": ["1_000000000000_000000000010"],
                        }
                    ],
                },
            )
            with self.assertRaisesRegex(GenerationManifestError, "noncanonical"):
                bucket_ids(path)

    def test_bucket_partitions_use_producer_chromosome_order(self) -> None:
        with TemporaryDirectory() as temporary:
            path = self.write(
                Path(temporary) / "buckets.json",
                {
                    "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
                    "buckets": [
                        {
                            "bucket_id": "00000",
                            "ordinal": 0,
                            "unit_ids": ["UNIT_A"],
                            "partition_ids": [
                                "1_000000000000_000000000010",
                                "10_000000000000_000000000010",
                                "2_000000000000_000000000010",
                            ],
                        }
                    ],
                },
            )
            self.assertEqual(("00000",), bucket_ids(path))

    def test_bucket_manifest_is_validated_once_per_published_version(self) -> None:
        with TemporaryDirectory() as temporary:
            path = self.write(
                Path(temporary) / "buckets.json",
                {
                    "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
                    "buckets": [
                        {
                            "bucket_id": "00000",
                            "ordinal": 0,
                            "unit_ids": ["UNIT_A"],
                            "partition_ids": ["1_000000000000_000000000010"],
                        }
                    ],
                },
            )
            with patch.object(
                generation_manifests,
                "_read_bucket_manifest_uncached",
                wraps=generation_manifests._read_bucket_manifest_uncached,
            ) as reader:
                self.assertEqual(("00000",), bucket_ids(path))
                self.assertEqual(
                    "00000", bucket_from_manifest(path, "00000").identifier
                )
                self.assertEqual(("00000",), bucket_ids(path))
                reader.assert_called_once()

                previous = path.stat()
                self.write(
                    path,
                    {
                        "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
                        "buckets": [
                            {
                                "bucket_id": "00000",
                                "ordinal": 0,
                                "unit_ids": ["UNIT_B"],
                                "partition_ids": ["1_000000000000_000000000010"],
                            }
                        ],
                    },
                )
                os.utime(
                    path,
                    ns=(previous.st_atime_ns, previous.st_mtime_ns + 1_000_000_000),
                )
                self.assertEqual(
                    ("UNIT_B",), bucket_from_manifest(path, "00000").unit_ids
                )
                self.assertEqual(2, reader.call_count)

    def test_partition_manifest_is_validated_once_per_published_version(self) -> None:
        with TemporaryDirectory() as temporary:
            path = self.write(
                Path(temporary) / "partitions.json",
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
            with patch.object(
                generation_manifests,
                "_read_partition_manifest_uncached",
                wraps=generation_manifests._read_partition_manifest_uncached,
            ) as reader:
                first = read_partition_manifest(path)
                second = read_partition_manifest(str(path))
                self.assertIs(first, second)
                reader.assert_called_once()

    def test_sparse_bucket_map_omits_nonoverlapping_partition_pairs(self) -> None:
        partitions = [
            GenomePartition("1_000000000000_000000000010", "1", 0, 10, 0),
            GenomePartition("1_000000000010_000000000020", "1", 10, 20, 1),
            GenomePartition("1_000000000020_000000000030", "1", 20, 30, 2),
        ]
        self.assertEqual(
            {
                0: (
                    "1_000000000000_000000000010",
                    "1_000000000010_000000000020",
                ),
                1: ("1_000000000020_000000000030",),
            },
            sparse_bucket_partitions(
                partitions,
                [("UNIT_A", "1", 8, 12), ("UNIT_B", "1", 25, 28)],
                {"UNIT_A": 0, "UNIT_B": 1},
                10,
            ),
        )

    def test_partitions_cover_only_windows_with_eligible_testing_units(self) -> None:
        with TemporaryDirectory() as temporary:
            fasta_index = Path(temporary) / "reference.fa.fai"
            fasta_index.write_text(
                "1\t30\t3\t30\t31\nX\t30\t37\t30\t31\n", encoding="utf-8"
            )

            partitions = build_partitions(
                fasta_index,
                [("1", 2, 8), ("1", 9, 12), ("X", 20, 21)],
                10,
            )

            self.assertEqual(
                [
                    ("1_000000000000_000000000010", "1", 0, 10),
                    ("1_000000000010_000000000020", "1", 10, 20),
                    ("X_000000000020_000000000030", "X", 20, 30),
                ],
                [
                    (
                        partition.identifier,
                        partition.chromosome,
                        partition.start,
                        partition.end,
                    )
                    for partition in partitions
                ],
            )

    def test_indel_halo_schedules_only_the_boundary_owner_partition(self) -> None:
        self.assertEqual(
            [("A", "1", 0, 3), ("B", "1", 0, 3)],
            indel_candidate_owner_intervals([("A", "1", 1, 3), ("B", "1", 1, 3)], 3),
        )
        partitions = [
            GenomePartition("1_000000000000_000000000010", "1", 0, 10, 0),
            GenomePartition("1_000000000010_000000000020", "1", 10, 20, 1),
            GenomePartition("1_000000000020_000000000030", "1", 20, 30, 2),
        ]
        callable_intervals = [("UNIT_A", "1", 10, 12)]
        owners = indel_candidate_owner_intervals(
            callable_intervals, max_reference_length=3
        )
        self.assertEqual(owners, [("UNIT_A", "1", 8, 12)])
        self.assertEqual(
            {
                0: (
                    "1_000000000000_000000000010",
                    "1_000000000010_000000000020",
                )
            },
            sparse_bucket_partitions(
                partitions,
                callable_intervals,
                {"UNIT_A": 0},
                10,
                candidate_owner_intervals=owners,
            ),
        )

    def test_partition_builder_includes_indel_owner_halo(self) -> None:
        with TemporaryDirectory() as temporary:
            fasta_index = Path(temporary) / "reference.fa.fai"
            fasta_index.write_text("1\t30\t3\t30\t31\n", encoding="utf-8")
            partitions = build_partitions(
                fasta_index,
                [("1", 10, 12)],
                10,
                candidate_owner_intervals=[
                    (chromosome, start, end)
                    for _, chromosome, start, end in indel_candidate_owner_intervals(
                        [("UNIT_A", "1", 10, 12)], max_reference_length=3
                    )
                ],
            )
        self.assertEqual(
            ["1_000000000000_000000000010", "1_000000000010_000000000020"],
            [partition.identifier for partition in partitions],
        )


if __name__ == "__main__":
    unittest.main()
