"""End-to-end N-by-M generation and ACAT reporting coverage."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from workflow_test_utils import SNAKEMAKE, command_output, run_workflow


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
class FullWorkflowFixtureTests(unittest.TestCase):
    def test_two_scored_partitions_write_raw_and_mid_p_acat_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory(prefix="portable-psap-acat-workflow-") as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow("tests/null_distribution/full_acat.run.yaml", outdir)
            self.assertEqual(0, result.returncode, command_output(result)[-6000:])
            cohort = outdir / "apply" / "acat" / "cohort"
            with gzip.open(cohort / "cohort_summary.tsv.gz", "rt", encoding="utf-8") as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(1, len(summary))
            self.assertEqual("2", summary[0]["contributing_partition_count"])
            self.assertEqual("upper_endpoint_unresolved", summary[0]["aggregation_status"])
            with gzip.open(cohort / "aggregation_diagnostics.tsv.gz", "rt", encoding="utf-8") as handle:
                diagnostics = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual("mid_p_exact_one_adjusted", diagnostics[0]["mid_p_aggregation_status"])
            self.assertGreater(float(diagnostics[0]["mid_p_ACAT_PSAP"]), 0)

    def test_two_partition_full_fixture_certifies_and_reports_all_partitions(self) -> None:
        with tempfile.TemporaryDirectory(prefix="portable-psap-full-workflow-") as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow("tests/null_distribution/full_partitioned.run.yaml", outdir)
            self.assertEqual(0, result.returncode, command_output(result)[-6000:])
            bundle = outdir / "generate" / "fixture" / "bundle"
            cohort = outdir / "apply" / "fixture" / "cohort"
            metadata = json.loads((bundle / "metadata.json").read_text())
            self.assertEqual("partitioned_null_bundle_v1", metadata["contract"])
            self.assertEqual(["cds", "intronic"], [entry["id"] for entry in metadata["partitions"]])
            self.assertEqual(10, len([key for key in metadata["components"] if key.startswith("lookup.")]))
            for partition in ("cds", "intronic"):
                self.assertTrue((bundle / f"resolved.{partition}.bed.gz").is_file())
                self.assertTrue((bundle / f"eligible.{partition}.bed.gz").is_file())
                self.assertTrue((cohort / "partitions" / f"{partition}.summary.tsv.gz").is_file())
            for name in ("cohort_summary.tsv.gz", "endpoint_diagnostics.tsv.gz", "aggregation_diagnostics.tsv.gz"):
                self.assertTrue((cohort / name).is_file())
            with gzip.open(cohort / "cohort_summary.tsv.gz", "rt", encoding="utf-8") as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(summary)
            self.assertEqual("ACAT", summary[0]["aggregation_method"])
            self.assertEqual("2", summary[0]["configured_partition_count"])
            with gzip.open(cohort / "endpoint_diagnostics.tsv.gz", "rt", encoding="utf-8") as handle:
                diagnostics = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(diagnostics)
            self.assertIn("exact_one_carrier_rate", diagnostics[0])
            resumed = run_workflow("tests/null_distribution/full_partitioned.run.yaml", outdir, dry_run=True)
            self.assertEqual(0, resumed.returncode, command_output(resumed)[-6000:])
            self.assertIn("Nothing to be done", command_output(resumed))


if __name__ == "__main__":
    unittest.main()
