"""Keep pre-migration main outputs as immutable historical evidence only."""

from __future__ import annotations

import json
from pathlib import Path
import csv
import gzip
import tempfile
import unittest

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


REFERENCE = ROOT / "tests" / "reference" / "partitioned_nulls" / "main_baseline"


class ReferenceWorkflowTests(unittest.TestCase):
    def test_pre_partition_baseline_is_documented_but_not_a_runtime_fixture(self) -> None:
        manifest = json.loads((REFERENCE / "manifest.json").read_text())
        self.assertTrue(manifest["cases"])
        for case in manifest["cases"]:
            with self.subTest(case=case["name"]):
                self.assertTrue((REFERENCE / case["name"] / "run.yaml").is_file())
                self.assertTrue(case["public_sha256"])

    @unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
    def test_n1_list_form_matches_frozen_role_rows_and_sample_psap(self) -> None:
        """The current N=1 partition path preserves pre-migration numerics."""
        with tempfile.TemporaryDirectory(prefix="portable-psap-n1-reference-") as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow("tests/null_distribution/full.run.yaml", outdir)
            self.assertEqual(0, result.returncode, command_output(result)[-6000:])
            bundle = outdir / "generate" / "fixture" / "bundle"

            def rows(path: Path) -> list[dict[str, str]]:
                with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                    return list(csv.DictReader(handle, delimiter="\t"))

            frozen = REFERENCE / "snv_hem"
            for current_name, frozen_name, unit in (
                ("lookup.cds.DOM-het.tsv.gz", "lookup_het.tsv.gz", "UNIT_OVERLAP"),
                ("lookup.cds.REC-hom.tsv.gz", "lookup_hom.tsv.gz", "UNIT_OVERLAP"),
                ("lookup.cds.X-linked-het.tsv.gz", "lookup_het.tsv.gz", "UNIT_X_MIXED_PAR"),
                ("lookup.cds.X-linked-hom.tsv.gz", "lookup_hom.tsv.gz", "UNIT_X_MIXED_PAR"),
                ("lookup.cds.X-linked-hem.tsv.gz", "lookup_hem_chrX.tsv.gz", "UNIT_X_MIXED_PAR"),
            ):
                with self.subTest(lookup=current_name):
                    current = rows(bundle / current_name)
                    previous = {
                        row["testing_unit_id"]: row
                        for row in rows(frozen / frozen_name)
                    }
                    self.assertEqual([unit], [row["testing_unit_id"] for row in current])
                    self.assertEqual(
                        len(current[0]), len(previous[unit]),
                        "the fixed-width score grid changed",
                    )
                    for column, value in current[0].items():
                        if column != "testing_unit_id":
                            self.assertAlmostEqual(float(value), float(previous[unit][column]), places=12)

            with gzip.open(
                outdir / "apply" / "fixture" / "scored" / "affected.gene.tsv.gz",
                "rt", encoding="utf-8", newline="",
            ) as handle:
                current_psaps = list(csv.DictReader(handle, delimiter="\t"))
            with (frozen / "sample_psaps.tsv").open(encoding="utf-8", newline="") as handle:
                frozen_psaps = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(1, len(current_psaps))
            self.assertEqual(1, len(frozen_psaps))
            self.assertEqual(frozen_psaps[0]["model"], current_psaps[0]["Dz.Model"])
            self.assertAlmostEqual(float(frozen_psaps[0]["PSAP"]), float(current_psaps[0]["PSAP"]), places=12)


if __name__ == "__main__":
    unittest.main()
