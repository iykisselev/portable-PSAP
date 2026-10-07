"""Fresh-output DAG baselines for every supported workflow mode."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import yaml

from workflow_test_utils import (
    ROOT,
    SNAKEMAKE,
    command_output,
    initial_job_count,
    run_workflow,
)


BASELINE_PATH = ROOT / "tests" / "workflow_baseline.json"


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow DAG tests")
class WorkflowDagBaselineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with BASELINE_PATH.open(encoding="utf-8") as handle:
            cls.baseline = json.load(handle)

    def test_fresh_mode_dags_match_the_recorded_baseline(self) -> None:
        for mode, expected in self.baseline.items():
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                prefix=f"portable-psap-{mode}-dag-"
            ) as temporary_dir:
                outdir = Path(temporary_dir) / "output"
                result = run_workflow(expected["manifest"], outdir, dry_run=True)
                output = command_output(result)
                self.assertEqual(0, result.returncode, output[-6000:])
                self.assertEqual(expected["initial_dry_run_jobs"], initial_job_count(result))

                resolved_output = ""
                checkpoint_target = expected.get("checkpoint_target")
                if checkpoint_target:
                    scope = run_workflow(
                        expected["manifest"], outdir, targets=(checkpoint_target,)
                    )
                    scope_output = command_output(scope)
                    self.assertEqual(0, scope.returncode, scope_output[-6000:])
                    resolved = run_workflow(
                        expected["manifest"], outdir, dry_run=True
                    )
                    resolved_output = command_output(resolved)
                    self.assertEqual(0, resolved.returncode, resolved_output[-6000:])
                    self.assertEqual(
                        expected["resolved_dry_run_jobs"], initial_job_count(resolved)
                    )

                for relative_path in (
                    expected["terminal_outputs"] + expected["important_outputs"]
                ):
                    self.assertIn(str(outdir / relative_path), output)
                for rule_name in expected["required_rules"]:
                    self.assertTrue(
                        any(
                            label in output + resolved_output
                            for label in (
                                f"rule {rule_name}:",
                                f"checkpoint {rule_name}:",
                            )
                        ),
                        msg=f"{mode} dry-run did not schedule {rule_name}",
                    )

    def test_removed_result_aliases_are_rejected(self) -> None:
        aliases = (
            ("tests/null_distribution/generate.run.yaml", "generate_results"),
            ("tests/null_distribution/generate.run.yaml", "null_distributions"),
            ("tests/smoke/run.yaml", "apply_results"),
            ("tests/null_distribution/full.run.yaml", "full_results"),
        )
        for manifest, target in aliases:
            with self.subTest(target=target), tempfile.TemporaryDirectory(
                prefix=f"portable-psap-{target}-"
            ) as temporary_dir:
                result = run_workflow(
                    manifest,
                    Path(temporary_dir) / "output",
                    dry_run=True,
                    targets=(target,),
                )
                output = command_output(result)
                self.assertNotEqual(0, result.returncode, output[-6000:])
                self.assertIn(f"No rule to produce {target}", output)

    def test_missing_enabled_indel_resource_reaches_validation_instead_of_dag_failure(self) -> None:
        with tempfile.TemporaryDirectory(prefix="portable-psap-missing-indel-") as temporary_dir:
            root = Path(temporary_dir)
            manifest = yaml.safe_load((ROOT / "tests" / "smoke" / "run.yaml").read_text())
            manifest["bindings"]["reference"]["cadd"]["indel"] = {
                "data": None,
                "index": None,
            }
            manifest_path = root / "missing-indel.yaml"
            manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
            result = run_workflow(manifest_path, root / "output", dry_run=True)
            output = command_output(result)
            self.assertEqual(0, result.returncode, output[-6000:])
            self.assertNotIn("missing input files: ''", output.lower())
            self.assertIn("rule validate_resources:", output)


if __name__ == "__main__":
    unittest.main()
