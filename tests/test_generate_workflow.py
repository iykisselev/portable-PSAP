"""End-to-end execution coverage for the mutation-rate generate fixture."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from workflow_test_utils import SNAKEMAKE, command_output, run_workflow


MUTATION_MANIFEST = "tests/null_distribution/mutation.run.yaml"
INDEL_MANIFEST = "tests/null_distribution/indel.run.yaml"


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow E2E tests")
class GenerateWorkflowFixtureTests(unittest.TestCase):
    def test_mutation_rate_fixture_publishes_a_certified_bundle(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-mutation-workflow-"
        ) as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow(MUTATION_MANIFEST, outdir)
            output = command_output(result)
            self.assertEqual(0, result.returncode, output[-6000:])

            bundle = outdir / "generate" / "mutation-fixture" / "bundle"
            metadata = json.loads(
                (bundle / "metadata.json").read_text(encoding="utf-8")
            )
            validation = json.loads(
                (bundle / "validation.json").read_text(encoding="utf-8")
            )
            self.assertEqual("generated_bundle_v2", metadata["contract"])
            self.assertEqual(
                "mutation_rate", metadata["model_assumptions"]["unobserved_policy"]
            )
            self.assertEqual(
                "de_novo_proxy", metadata["model_assumptions"]["mutation_conversion"]
            )
            self.assertEqual(
                "sex_averaged_proxy",
                metadata["model_assumptions"]["x_chromosome"]["mutation_rate_policy"],
            )
            self.assertIn("population_af_zero_imputations", metadata["components"])
            self.assertTrue(
                (bundle / "population_af_zero_imputations.tsv.gz").is_file()
            )
            self.assertEqual("ok", validation["status"])
            self.assertTrue((bundle / "validation.ok").is_file())

    def test_indel_fixture_publishes_a_certified_v3_bundle(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-indel-workflow-"
        ) as temporary:
            outdir = Path(temporary) / "output"
            result = run_workflow(INDEL_MANIFEST, outdir)
            output = command_output(result)
            self.assertEqual(0, result.returncode, output[-6000:])

            bundle = outdir / "generate" / "indel-fixture" / "bundle"
            metadata = json.loads(
                (bundle / "metadata.json").read_text(encoding="utf-8")
            )
            validation = json.loads(
                (bundle / "validation.json").read_text(encoding="utf-8")
            )
            self.assertEqual("generated_bundle_v3", metadata["contract"])
            self.assertEqual(
                ["snv", "indel"], metadata["model_assumptions"]["variant_classes"]
            )
            self.assertIn("score_snv", metadata["resource_identities"])
            self.assertIn("score_indel", metadata["resource_identities"])
            self.assertIn("candidate_audit", metadata["components"])
            self.assertEqual({"snv", "indel"}, set(metadata["variant_class_counts"]))
            self.assertTrue((bundle / "candidate_audit.tsv.gz").is_file())
            self.assertEqual("ok", validation["status"])


if __name__ == "__main__":
    unittest.main()
