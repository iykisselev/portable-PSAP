"""Ensure every documented run-mode demonstration still builds a DAG."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


EXAMPLES = {
    "apply": "runs/apply.example.yaml",
    "generate": "runs/generate.example.yaml",
    "full": "runs/full.example.yaml",
}


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for example-manifest tests")
class RunExampleTests(unittest.TestCase):
    def test_every_documented_example_constructs_its_mode_dag(self) -> None:
        for mode, manifest in EXAMPLES.items():
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                prefix=f"portable-psap-{mode}-example-"
            ) as temporary:
                result = run_workflow(manifest, Path(temporary) / "output", dry_run=True)
                output = command_output(result)
                self.assertEqual(0, result.returncode, output[-6000:])


if __name__ == "__main__":
    unittest.main()
