"""End-to-end coverage for the external-lookup apply workflow fixture."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import yaml

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


def expected_digests() -> dict[str, str]:
    """Read the stable digest manifest for the external-lookup reports."""
    digests = {}
    with (ROOT / "tests" / "smoke" / "expected_digests.sha256").open(
        encoding="utf-8"
    ) as handle:
        for line in handle:
            digest, filename = line.split()
            digests[filename] = digest
    return digests


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
class ApplyWorkflowFixtureTests(unittest.TestCase):
    def test_external_lookup_apply_fixture_completes_with_stable_reports(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-apply-workflow-"
        ) as temporary_dir:
            outdir = Path(temporary_dir) / "output"
            manifest = yaml.safe_load((ROOT / "tests/smoke/run.yaml").read_text())
            # Repeated exclusions must be subtracted once, never treated as an
            # inclusion mask by the renamed primary input.
            manifest["bindings"]["reference"]["exclusion"] = dict(
                manifest["bindings"]["reference"]["exclusions"][
                    "encode-hg38-blacklist-v2"
                ]
            )
            manifest_path = Path(temporary_dir) / "run.yaml"
            manifest_path.write_text(yaml.safe_dump(manifest))
            result = run_workflow(manifest_path, outdir)
            output = command_output(result)
            self.assertEqual(0, result.returncode, output[-6000:])

            apply_dir = outdir / "apply" / "smoke"
            cohort = apply_dir / "cohort"
            self.assertFalse((apply_dir / "standardized").exists())
            self.assertFalse((apply_dir / "sex_classified").exists())
            self.assertFalse(
                any((apply_dir / "raw").rglob("*.tsv.gz")),
                "temporary extracted genotype tables were retained",
            )
            self.assertFalse(
                any((apply_dir / "raw").rglob("*.vcf.gz")),
                "temporary split genotype VCFs were retained",
            )
            for role in ("het", "hom", "hem_chrX"):
                self.assertTrue(
                    (apply_dir / "prepared" / "lookups" / f"{role}.rds").is_file(),
                    f"missing prepared lookup for {role}",
                )
            for filename, digest in expected_digests().items():
                report = cohort / filename
                self.assertTrue(report.is_file(), f"missing workflow output: {report}")
                self.assertEqual(
                    digest, sha256(report), f"changed fixture digest: {filename}"
                )

            with (apply_dir / "validation" / "model.json").open(
                encoding="utf-8"
            ) as handle:
                model_validation = json.load(handle)
            self.assertEqual("ok", model_validation["status"])
            self.assertTrue(model_validation["summary"]["hemizygote_model"])
            self.assertEqual(
                "external_legacy_v0", model_validation["summary"]["lookup_contract"]
            )

            with gzip.open(
                cohort / "cohort_candidates.tsv.gz", "rt", encoding="utf-8"
            ) as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(candidates, "apply fixture produced no candidates")
            self.assertIn("X-linked-hem", {row["model"] for row in candidates})


if __name__ == "__main__":
    unittest.main()
