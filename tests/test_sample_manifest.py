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
    read_sample_manifest,
    read_sample_manifest_index,
    sample_ids_by_phenotype,
    sample_metadata,
    sample_sexes_by_phenotype,
)


class SampleManifestTests(unittest.TestCase):
    def test_role_queries_preserve_manifest_order(self) -> None:
        manifest = read_sample_manifest(ROOT / "tests" / "smoke" / "data" / "samples.tsv")

        self.assertEqual(["affected"], sample_ids_by_phenotype(manifest, "2"))
        self.assertEqual(["unaffected"], sample_ids_by_phenotype(manifest, "1"))
        self.assertEqual(["2"], sample_sexes_by_phenotype(manifest, "1"))
        self.assertEqual("1", sample_metadata(manifest, "affected")["sex"])

    def test_expected_role_is_checked_against_the_loaded_manifest(self) -> None:
        manifest = read_sample_manifest(ROOT / "tests" / "smoke" / "data" / "samples.tsv")

        with self.assertRaisesRegex(ValueError, "expected 'affected'"):
            sample_metadata(manifest, "unaffected", expected_phenotype="affected")

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
