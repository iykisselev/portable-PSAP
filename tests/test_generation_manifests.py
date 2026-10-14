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
    WINDOW_MANIFEST_SCHEMA_VERSION,
    GenerationManifestError,
    WindowManifestEntry,
    bucket_from_manifest,
    bucket_ids,
    window_ids,
    read_window_manifest,
)
from generation_windows import (  # noqa: E402
    GenomeWindow,
    build_windows,
    indel_candidate_owner_intervals,
    read_genome_windows,
)
from generation_scope import sparse_bucket_windows  # noqa: E402


class GenerationManifestTests(unittest.TestCase):
    def write(self, path: Path, payload: dict) -> Path:
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_pre_rename_manifests_require_regeneration(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            windows = self.write(
                directory / "old-windows.json",
                {"schema_version": 1, "partitions": []},
            )
            buckets = self.write(
                directory / "old-buckets.json",
                {"schema_version": 2, "buckets": []},
            )
            with self.assertRaisesRegex(GenerationManifestError, "Unsupported"):
                read_window_manifest(windows)
            with self.assertRaisesRegex(GenerationManifestError, "Unsupported"):
                bucket_ids(buckets)

    def test_old_window_fields_are_not_accepted_under_the_new_version(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            windows = self.write(
                directory / "windows.json",
                {"schema_version": WINDOW_MANIFEST_SCHEMA_VERSION, "partitions": []},
            )
            buckets = self.write(
                directory / "buckets.json",
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
            with self.assertRaisesRegex(GenerationManifestError, "windows array"):
                read_window_manifest(windows)
            with self.assertRaisesRegex(GenerationManifestError, "invalid entry"):
                bucket_ids(buckets)

    def test_validated_manifests_publish_canonical_dynamic_wildcards(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            windows = self.write(
                directory / "windows.json",
                {
                    "schema_version": WINDOW_MANIFEST_SCHEMA_VERSION,
                    "windows": [
                        {
                            "window_id": "1_000000000000_000000000010",
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
                            "window_ids": ["1_000000000000_000000000010"],
                        }
                    ],
                },
            )

            self.assertEqual(("1_000000000000_000000000010",), window_ids(windows))
            self.assertIsInstance(read_window_manifest(windows)[0], WindowManifestEntry)
            self.assertIsInstance(read_genome_windows(windows)[0], GenomeWindow)
            self.assertEqual(("00000",), bucket_ids(buckets))
            self.assertEqual(
                ("UNIT_A", "UNIT_B"), bucket_from_manifest(buckets, "00000").unit_ids
            )
            self.assertEqual(
                ("1_000000000000_000000000010",),
                bucket_from_manifest(buckets, "00000").window_ids,
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
                            "window_ids": ["1_000000000000_000000000010"],
                        }
                    ],
                },
            )
            with self.assertRaisesRegex(GenerationManifestError, "noncanonical"):
                bucket_ids(path)

    def test_bucket_windows_use_producer_chromosome_order(self) -> None:
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
                            "window_ids": [
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
                            "window_ids": ["1_000000000000_000000000010"],
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
                                "window_ids": ["1_000000000000_000000000010"],
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

    def test_window_manifest_is_validated_once_per_published_version(self) -> None:
        with TemporaryDirectory() as temporary:
            path = self.write(
                Path(temporary) / "windows.json",
                {
                    "schema_version": WINDOW_MANIFEST_SCHEMA_VERSION,
                    "windows": [
                        {
                            "window_id": "1_000000000000_000000000010",
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
                "_read_window_manifest_uncached",
                wraps=generation_manifests._read_window_manifest_uncached,
            ) as reader:
                first = read_window_manifest(path)
                second = read_window_manifest(str(path))
                self.assertIs(first, second)
                reader.assert_called_once()

    def test_sparse_bucket_map_omits_nonoverlapping_window_pairs(self) -> None:
        windows = [
            GenomeWindow("1_000000000000_000000000010", "1", 0, 10, 0),
            GenomeWindow("1_000000000010_000000000020", "1", 10, 20, 1),
            GenomeWindow("1_000000000020_000000000030", "1", 20, 30, 2),
        ]
        self.assertEqual(
            {
                0: (
                    "1_000000000000_000000000010",
                    "1_000000000010_000000000020",
                ),
                1: ("1_000000000020_000000000030",),
            },
            sparse_bucket_windows(
                windows,
                [("UNIT_A", "1", 8, 12), ("UNIT_B", "1", 25, 28)],
                {"UNIT_A": 0, "UNIT_B": 1},
                10,
            ),
        )

    def test_windows_cover_only_windows_with_eligible_testing_units(self) -> None:
        with TemporaryDirectory() as temporary:
            fasta_index = Path(temporary) / "reference.fa.fai"
            fasta_index.write_text(
                "1\t30\t3\t30\t31\nX\t30\t37\t30\t31\n", encoding="utf-8"
            )

            windows = build_windows(
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
                        window.identifier,
                        window.chromosome,
                        window.start,
                        window.end,
                    )
                    for window in windows
                ],
            )

    def test_indel_halo_schedules_only_the_boundary_owner_window(self) -> None:
        self.assertEqual(
            [("A", "1", 0, 3), ("B", "1", 0, 3)],
            indel_candidate_owner_intervals([("A", "1", 1, 3), ("B", "1", 1, 3)], 3),
        )
        windows = [
            GenomeWindow("1_000000000000_000000000010", "1", 0, 10, 0),
            GenomeWindow("1_000000000010_000000000020", "1", 10, 20, 1),
            GenomeWindow("1_000000000020_000000000030", "1", 20, 30, 2),
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
            sparse_bucket_windows(
                windows,
                callable_intervals,
                {"UNIT_A": 0},
                10,
                candidate_owner_intervals=owners,
            ),
        )

    def test_window_builder_includes_indel_owner_halo(self) -> None:
        with TemporaryDirectory() as temporary:
            fasta_index = Path(temporary) / "reference.fa.fai"
            fasta_index.write_text("1\t30\t3\t30\t31\n", encoding="utf-8")
            windows = build_windows(
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
            [window.identifier for window in windows],
        )


if __name__ == "__main__":
    unittest.main()
