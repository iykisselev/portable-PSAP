"""Adversarial apply-mode checks for certified generated-bundle handoff."""

from __future__ import annotations

import copy
from pathlib import Path
import shutil
import tempfile
import unittest

import yaml

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


GENERATE_MANIFEST = "tests/null_distribution/generate.run.yaml"


def bundle_lookup_source(bundle: Path) -> dict[str, object]:
    """Build the only supported apply-mode binding shape for one bundle."""
    return {
        "type": "bundle",
        "contract": "generated_bundle_v2",
        "tables": {
            "het": str(bundle / "lookup_het.tsv.gz"),
            "hom": str(bundle / "lookup_hom.tsv.gz"),
            "hem_chrX": str(bundle / "lookup_hem_chrX.tsv.gz"),
        },
        "metadata": str(bundle / "metadata.json"),
        "validation": str(bundle / "validation.json"),
        "marker": str(bundle / "validation.ok"),
        "eligible_units": str(bundle / "eligible_testing_units.bed.gz"),
        "calibration_mask": str(bundle / "calibration_regions.bed.gz"),
        "excluded_alleles": str(bundle / "excluded_alleles.tsv.gz"),
        "excluded_units": str(bundle / "excluded_units.tsv.gz"),
    }


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow resilience tests")
class BundleFailureWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory(
            prefix="portable-psap-bundle-failure-"
        )
        cls.root = Path(cls.temporary.name)
        result = run_workflow(GENERATE_MANIFEST, cls.root / "generated")
        output = command_output(result)
        if result.returncode:
            raise AssertionError(output[-6000:])
        cls.bundle = cls.root / "generated" / "generate" / "fixture" / "bundle"
        cls.base_manifest = yaml.safe_load(
            (ROOT / "tests" / "smoke" / "run.yaml").read_text(encoding="utf-8")
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def apply_manifest(self, bundle: Path, name: str) -> Path:
        """Write an isolated apply manifest that binds exactly one certificate."""
        manifest = copy.deepcopy(self.base_manifest)
        manifest["workflow"] = {"mode": "apply"}
        manifest["run"] = {"apply_prefix": name}
        manifest["bindings"]["reference"]["lookup_source"] = bundle_lookup_source(
            bundle
        )
        path = self.root / f"{name}.yaml"
        path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
        return path

    def test_corrupt_component_fails_certificate_validation_before_apply_work(
        self,
    ) -> None:
        bundle = self.root / "corrupt-bundle"
        shutil.copytree(self.bundle, bundle)
        with (bundle / "lookup_het.tsv.gz").open("ab") as handle:
            handle.write(b"corrupt")
        result = run_workflow(
            self.apply_manifest(bundle, "corrupt"), self.root / "corrupt-output"
        )
        output = command_output(result)
        self.assertNotEqual(0, result.returncode)
        self.assertRegex(
            output,
            r"Generated bundle component lookup_het has the wrong (size|checksum)",
        )
        self.assertNotIn("rule annotate_variants:", output)

    def test_missing_declared_component_fails_before_apply_work(self) -> None:
        bundle = self.root / "incomplete-bundle"
        shutil.copytree(self.bundle, bundle)
        (bundle / "excluded_units.tsv.gz").unlink()
        result = run_workflow(
            self.apply_manifest(bundle, "incomplete"), self.root / "incomplete-output"
        )
        output = command_output(result)
        self.assertNotEqual(0, result.returncode)
        self.assertIn("excluded_units.tsv.gz", output)
        self.assertNotIn("rule annotate_variants:", output)


if __name__ == "__main__":
    unittest.main()
