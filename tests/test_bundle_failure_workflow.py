"""Adversarial apply checks for certified partitioned-bundle handoff."""

from __future__ import annotations

import copy
from pathlib import Path
import shutil
import tempfile
import unittest

import yaml

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow resilience tests")
class BundleFailureWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory(prefix="portable-psap-bundle-failure-")
        cls.root = Path(cls.temporary.name)
        result = run_workflow("tests/null_distribution/generate.run.yaml", cls.root / "generated")
        if result.returncode:
            raise AssertionError(command_output(result)[-6000:])
        cls.bundle = cls.root / "generated" / "generate" / "fixture" / "bundle"
        cls.base_manifest = yaml.safe_load((ROOT / "tests/null_distribution/full.run.yaml").read_text())

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def apply_manifest(self, bundle: Path, name: str) -> Path:
        manifest = copy.deepcopy(self.base_manifest)
        manifest["workflow"] = {"mode": "apply"}
        manifest["run"] = {"apply_prefix": name}
        manifest.pop("generation")
        manifest["reference"].pop("testing_unit_partitions")
        reference = manifest["bindings"]["reference"]
        for key in ("exclusion", "exclusions", "testing_unit_partitions"):
            reference.pop(key)
        reference["lookup_source"] = {"type": "bundle", "directory": str(bundle)}
        path = self.root / f"{name}.yaml"
        path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
        return path

    def test_corrupt_lookup_fails_certificate_before_annotation(self) -> None:
        bundle = self.root / "corrupt-bundle"
        shutil.copytree(self.bundle, bundle)
        with (bundle / "lookup.cds.DOM-het.tsv.gz").open("ab") as handle:
            handle.write(b"corrupt")
        result = run_workflow(self.apply_manifest(bundle, "corrupt"), self.root / "corrupt-output")
        output = command_output(result)
        self.assertNotEqual(0, result.returncode)
        self.assertRegex(output, r"lookup\.cds\.DOM-het has the wrong (size|checksum)")
        self.assertNotIn("rule annotate_variants:", output)

    def test_missing_certified_partition_component_fails_before_annotation(self) -> None:
        bundle = self.root / "incomplete-bundle"
        shutil.copytree(self.bundle, bundle)
        (bundle / "eligible.cds.bed.gz").unlink()
        result = run_workflow(self.apply_manifest(bundle, "incomplete"), self.root / "incomplete-output")
        output = command_output(result)
        self.assertNotEqual(0, result.returncode)
        self.assertIn("eligible.cds.bed.gz", output)
        self.assertNotIn("rule annotate_variants:", output)


if __name__ == "__main__":
    unittest.main()
