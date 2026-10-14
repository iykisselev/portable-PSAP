"""Integrity and scientific coverage of the immutable pre-feature reference."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import unittest


REFERENCE = Path(__file__).parent / "reference" / "partitioned_nulls" / "main_baseline"


class PartitionedReferenceTests(unittest.TestCase):
    def test_frozen_files_match_the_capture_manifest(self):
        manifest = json.loads((REFERENCE / "manifest.json").read_text())
        self.assertEqual(40, len(manifest["commit"]))
        self.assertEqual(6, len(manifest["cases"]))
        for relative, expected in manifest["reference_sha256"].items():
            with self.subTest(path=relative):
                self.assertEqual(
                    expected,
                    hashlib.sha256((REFERENCE / relative).read_bytes()).hexdigest(),
                )

    def test_reference_covers_all_models_and_a_retained_indel(self):
        with (REFERENCE / "models_snv_hem" / "sample_psaps.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            {"DOM-het", "REC-hom", "X-linked-hom", "X-linked-het", "X-linked-hem"},
            {row["model"] for row in rows},
        )
        with (REFERENCE / "models_snv_no_hem" / "sample_psaps.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(4, len(rows))
        self.assertNotIn("X-linked-hem", {row["model"] for row in rows})
        with gzip.open(
            REFERENCE / "models_indel_hem" / "male.scored.tsv.gz", "rt"
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        indels = [row for row in rows if row["Key"] == "1:3-3_A>AT"]
        self.assertEqual(1, len(indels))
        self.assertEqual("DOM-het", indels[0]["Dz.Model"])
        self.assertGreater(float(indels[0]["PSAP"]), 0)
        self.assertLess(float(indels[0]["PSAP"]), 1)


if __name__ == "__main__":
    unittest.main()
