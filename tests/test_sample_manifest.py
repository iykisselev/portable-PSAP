"""Contracts for reusable apply sample-manifest helpers."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from sample_manifest import (  # noqa: E402
    configured_sample_index,
    read_sample_manifest_index,
)
from sample_metadata import SampleRow, normalize_sample  # noqa: E402


class SampleManifestTests(unittest.TestCase):
    def test_shared_normalizer_defines_safe_manifest_rows(self) -> None:
        sample = normalize_sample(SampleRow(7, " F ", " affected ", " 1 ", " 2 "))
        self.assertTrue(sample.is_safe)
        self.assertEqual(
            {"family_id": "F", "pid": "affected", "sex": "1", "phenotype": "affected"},
            sample.manifest_row(),
        )
        unsafe = normalize_sample(SampleRow(8, "F", "bad/name", "2", "1"))
        self.assertFalse(unsafe.is_safe)

    def test_configured_index_is_immutable_and_keeps_safe_manifest_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            samples = Path(temporary_dir) / "samples.tsv"
            samples.write_text(
                "family_id\tpid\tsex\tphenotype\n"
                "F1\taffected\t1\t2\n"
                "F1\tunknown\t2\t0\n"
                "F1\tunaffected\t2\t1\n"
                "F1\taffected\t2\t1\n"
                "F1\tunsafe/name\t1\t2\n",
                encoding="utf-8",
            )

            index = configured_sample_index(samples, "samples")

            self.assertEqual(["affected"], index.sample_ids("affected"))
            self.assertEqual(["unaffected"], index.sample_ids("unaffected"))
            self.assertEqual("1", index.sample("affected").sex)
            with self.assertRaises(TypeError):
                index.by_pid["new"] = index.sample("affected")

    def test_canonical_manifest_index_uses_the_published_rows(self) -> None:
        index = read_sample_manifest_index(
            ROOT / "tests" / "smoke" / "data" / "samples.tsv"
        )

        self.assertEqual(["affected"], index.sample_ids("2"))
        self.assertEqual(["unaffected"], index.sample_ids("1"))


if __name__ == "__main__":
    unittest.main()
