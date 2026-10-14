"""Compare the window refactor against real outputs captured from unchanged main."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


REFERENCE = ROOT / "tests" / "reference" / "partitioned_nulls" / "main_baseline"


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for reference workflow tests")
class ReferenceWorkflowTests(unittest.TestCase):
    def test_public_outputs_are_identical_to_the_captured_main_outputs(self):
        manifest = json.loads((REFERENCE / "manifest.json").read_text())
        for case in manifest["cases"]:
            with (
                self.subTest(case=case["name"]),
                tempfile.TemporaryDirectory(
                    prefix="psap-window-reference-"
                ) as temporary,
            ):
                output = Path(temporary) / "output"
                result = run_workflow(REFERENCE / case["name"] / "run.yaml", output)
                self.assertEqual(0, result.returncode, command_output(result)[-6000:])
                for relative, expected in case["public_sha256"].items():
                    with self.subTest(component=relative):
                        self.assertEqual(
                            expected,
                            hashlib.sha256(
                                (output / relative).read_bytes()
                            ).hexdigest(),
                        )


if __name__ == "__main__":
    unittest.main()
