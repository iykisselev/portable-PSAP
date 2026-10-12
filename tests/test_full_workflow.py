"""End-to-end execution exclusion for the staged full workflow fixture."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import yaml

from workflow_test_utils import (
    ROOT,
    SNAKEMAKE,
    command_output,
    completed_step_count,
    planned_job_count,
    run_workflow,
)

BASELINE_PATH = ROOT / "tests" / "workflow_baseline.json"
FULL_MANIFEST = "tests/null_distribution/full.run.yaml"
FULL_INDEL_MANIFEST = "tests/null_distribution/full_indel.run.yaml"


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
class FullWorkflowFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with BASELINE_PATH.open(encoding="utf-8") as handle:
            cls.baseline = json.load(handle)["full"]

    def test_full_fixture_generates_and_applies_hem_bundle(self) -> None:
        """Execute the complete DAG and verify its generated bundle is consumed."""
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-full-workflow-"
        ) as temporary_dir:
            temporary_root = Path(temporary_dir)
            outdir = temporary_root / "output"
            result = run_workflow(FULL_MANIFEST, outdir)
            if result.returncode:
                output = command_output(result)[-6000:]
                self.fail(
                    f"full fixture failed with exit code {result.returncode}:\n{output}"
                )
            self.assertEqual(
                self.baseline["fresh_execution_steps"], completed_step_count(result)
            )

            bundle = outdir / "generate" / "fixture" / "bundle"
            cohort = outdir / "apply" / "fixture" / "cohort"
            for path in (
                bundle / "lookup_het.tsv.gz",
                bundle / "lookup_hom.tsv.gz",
                bundle / "lookup_hem_chrX.tsv.gz",
                bundle / "metadata.json",
                bundle / "validation.json",
                bundle / "validation.ok",
                outdir / "apply" / "fixture" / "validation" / "model.json",
                cohort / "cohort_candidates.tsv.gz",
                cohort / "cohort_summary.tsv.gz",
                cohort / "control_variants.tsv.gz",
            ):
                self.assertTrue(path.is_file(), f"missing workflow output: {path}")
                self.assertGreater(
                    path.stat().st_size, 0, f"empty workflow output: {path}"
                )

            with (bundle / "metadata.json").open(encoding="utf-8") as handle:
                metadata = json.load(handle)
            self.assertEqual("generated_bundle_v2", metadata["contract"])
            self.assertIn("lookup_hem_chrX", metadata["components"])
            self.assertIn("hem_chrX", metadata["model_assumptions"]["roles"])

            with (ROOT / "tests" / "null_distribution" / "full.run.yaml").open(
                encoding="utf-8"
            ) as handle:
                shared_reference = yaml.safe_load(handle)["bindings"]["reference"]
            resource_paths = metadata["resource_identities"]
            self.assertEqual(
                shared_reference["fasta"]["data"], resource_paths["fasta"]["path"]
            )
            self.assertEqual(
                shared_reference["cadd"]["snv"]["data"],
                resource_paths["score_snv"]["path"],
            )
            self.assertEqual(
                shared_reference["testing_unit_bed"]["data"],
                resource_paths["testing_units"]["path"],
            )
            self.assertEqual(
                shared_reference["exclusion"]["data"],
                resource_paths["exclusion"]["path"],
            )
            self.assertEqual(
                shared_reference["exclusions"]["encode-hg38-blacklist-v2"]["data"],
                resource_paths["encode-hg38-blacklist-v2"]["path"],
            )
            self.assertEqual(
                shared_reference["exclusions"]["par-regions"]["data"],
                resource_paths["par_regions"]["path"],
            )

            with (outdir / "apply" / "fixture" / "validation" / "model.json").open(
                encoding="utf-8"
            ) as handle:
                model_validation = json.load(handle)
            self.assertEqual("ok", model_validation["status"])
            self.assertTrue(model_validation["summary"]["hemizygote_model"])
            self.assertEqual(
                "generated_bundle_v2", model_validation["summary"]["lookup_contract"]
            )
            self.assertIn(
                str(bundle / "lookup_hem_chrX.tsv.gz"),
                model_validation["summary"]["lookup_files_checked"],
            )

            with gzip.open(
                cohort / "cohort_summary.tsv.gz", "rt", encoding="utf-8"
            ) as handle:
                header = handle.readline().rstrip("\n").split("\t")
            self.assertIn("minimum_PSAP", header)

            resumed = run_workflow(FULL_MANIFEST, outdir, dry_run=True)
            resumed_output = command_output(resumed)
            self.assertEqual(0, resumed.returncode, resumed_output[-6000:])
            self.assertIn("Nothing to be done", resumed_output)
            self.assertEqual(
                self.baseline["resume_dry_run_jobs"], planned_job_count(resumed)
            )

    def test_full_indel_fixture_consumes_the_generated_v3_bundle(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-full-indel-workflow-"
        ) as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow(FULL_INDEL_MANIFEST, outdir)
            output = command_output(result)
            self.assertEqual(0, result.returncode, output[-6000:])

            bundle = outdir / "generate" / "indel-fixture" / "bundle"
            model = outdir / "apply" / "indel-fixture" / "validation" / "model.json"
            metadata = json.loads(
                (bundle / "metadata.json").read_text(encoding="utf-8")
            )
            model_validation = json.loads(model.read_text(encoding="utf-8"))
            self.assertEqual("generated_bundle_v3", metadata["contract"])
            self.assertEqual(
                "generated_bundle_v3", model_validation["summary"]["lookup_contract"]
            )
            self.assertTrue((bundle / "validation.ok").is_file())


if __name__ == "__main__":
    unittest.main()
