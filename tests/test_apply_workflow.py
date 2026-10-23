"""End-to-end coverage for standalone apply from a certified partitioned bundle."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import yaml

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
class ApplyWorkflowFixtureTests(unittest.TestCase):
    def test_bundle_apply_emits_partition_and_omnibus_reports(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-apply-workflow-"
        ) as temporary:
            root = Path(temporary)
            generated = root / "generated"
            generated_result = run_workflow(
                "tests/null_distribution/generate.run.yaml", generated
            )
            self.assertEqual(
                0, generated_result.returncode, command_output(generated_result)[-6000:]
            )
            bundle = generated / "generate" / "fixture" / "bundle"
            manifest = yaml.safe_load(
                (ROOT / "tests/null_distribution/full.run.yaml").read_text()
            )
            manifest["workflow"] = {"mode": "apply"}
            manifest["run"].pop("generation_prefix")
            manifest.pop("generation")
            manifest["reference"].pop("testing_unit_partitions")
            reference = manifest["bindings"]["reference"]
            for name in ("exclusion", "exclusions", "testing_unit_partitions"):
                reference.pop(name)
            reference["lookup_source"] = {"type": "bundle", "directory": str(bundle)}
            path = root / "apply.yaml"
            path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
            result = run_workflow(path, root / "applied")
            self.assertEqual(0, result.returncode, command_output(result)[-6000:])
            cohort = root / "applied" / "apply" / "fixture" / "cohort"
            for name in (
                "cohort_summary.tsv.gz",
                "endpoint_diagnostics.tsv.gz",
                "aggregation_diagnostics.tsv.gz",
                "partitions/cds.summary.tsv.gz",
            ):
                self.assertTrue((cohort / name).is_file(), name)
            with gzip.open(
                cohort / "cohort_summary.tsv.gz", "rt", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(rows)
            self.assertEqual("1", rows[0]["configured_partition_count"])
            self.assertEqual(
                float(rows[0]["PSAP"]),
                float(json.loads(rows[0]["partition_psaps_json"])[0]["psap"]),
            )
            metadata = json.loads((cohort / "report_metadata.json").read_text())
            self.assertEqual("raw_survival_v1", metadata["aggregation_input_policy"])


if __name__ == "__main__":
    unittest.main()
