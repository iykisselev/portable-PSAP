"""Ensure every documented run-mode demonstration still builds a DAG."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import yaml

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


EXAMPLES = {
    "apply": "runs/apply.example.yaml",
    "generate": "runs/generate.example.yaml",
    "full": "runs/full.example.yaml",
}


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for example-manifest tests")
class RunExampleTests(unittest.TestCase):
    def test_combined_examples_bind_separate_indel_cadd_resources(self) -> None:
        generate = yaml.safe_load((ROOT / "runs/generate.example.yaml").read_text())
        full = yaml.safe_load((ROOT / "runs/full.example.yaml").read_text())

        self.assertTrue(generate["model"]["include_indels"])
        self.assertEqual(
            ["snv", "indel"], generate["generation_settings"]["variant_classes"]
        )
        self.assertIn("score_indel", generate["generation"])

        self.assertTrue(full["model"]["include_indels"])
        self.assertEqual(
            ["snv", "indel"], full["generation_settings"]["variant_classes"]
        )
        self.assertIn("indel", full["bindings"]["reference"]["cadd"])

    def test_unindexed_production_generation_fails_during_dag_construction(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-index-required-"
        ) as temporary:
            directory = Path(temporary)
            manifest = yaml.safe_load((ROOT / "runs/generate.example.yaml").read_text())
            manifest["generation"]["execution"].pop("allow_unindexed_fixture_inputs")
            manifest_path = directory / "run.yaml"
            manifest_path.write_text(yaml.safe_dump(manifest))
            result = run_workflow(manifest_path, directory / "output", dry_run=True)
            output = command_output(result)
            self.assertNotEqual(0, result.returncode)
            self.assertIn("generation.population_vcf.index is required", output)
            self.assertNotIn("Job stats:", output)

    def test_every_documented_example_constructs_its_mode_dag(self) -> None:
        for mode, manifest in EXAMPLES.items():
            with (
                self.subTest(mode=mode),
                tempfile.TemporaryDirectory(
                    prefix=f"portable-psap-{mode}-example-"
                ) as temporary,
            ):
                result = run_workflow(
                    manifest, Path(temporary) / "output", dry_run=True
                )
                output = command_output(result)
                self.assertEqual(0, result.returncode, output[-6000:])


if __name__ == "__main__":
    unittest.main()
