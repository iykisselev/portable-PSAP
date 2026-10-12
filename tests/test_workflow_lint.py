"""Ensure every workflow mode has only reviewed Snakemake lint findings."""

from __future__ import annotations

import json
from pathlib import Path
import re
import tempfile
import unittest

from workflow_test_utils import ROOT, SNAKEMAKE, command_output, run_workflow


ALLOWLIST_PATH = ROOT / "tests" / "workflow_lint_allowlist.json"
BASELINE_PATH = ROOT / "tests" / "workflow_baseline.json"
LINT_HEADER = re.compile(r"^Lints for snakefile (.+):$")
LINT_MESSAGE = re.compile(r"^    \* (.+)$")


def lint_messages(output: str) -> dict[str, list[str]]:
    """Extract lint headings by repository-relative snakefile path."""
    findings: dict[str, list[str]] = {}
    lint_path: str | None = None
    for line in output.splitlines():
        header = LINT_HEADER.match(line)
        if header:
            lint_path = Path(header.group(1)).resolve().relative_to(ROOT).as_posix()
            findings[lint_path] = []
            continue

        message = LINT_MESSAGE.match(line)
        if message and lint_path is not None:
            findings[lint_path].append(message.group(1))
    return findings


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for workflow lint tests")
class WorkflowLintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with BASELINE_PATH.open(encoding="utf-8") as handle:
            cls.baseline = json.load(handle)
        with ALLOWLIST_PATH.open(encoding="utf-8") as handle:
            cls.allowlist = json.load(handle)

    def test_lint_findings_are_explicitly_allowlisted_for_each_mode(self) -> None:
        for mode, expected in self.baseline.items():
            with (
                self.subTest(mode=mode),
                tempfile.TemporaryDirectory(
                    prefix=f"portable-psap-{mode}-lint-"
                ) as temporary_dir,
            ):
                result = run_workflow(
                    expected["manifest"], Path(temporary_dir) / "output", lint=True
                )
                output = command_output(result)
                observed = lint_messages(output)
                expected_findings = {
                    path: details["messages"]
                    for path, details in self.allowlist[mode].items()
                }
                self.assertEqual(expected_findings, observed, output[-6000:])
                self.assertEqual(
                    bool(expected_findings), bool(result.returncode), output[-6000:]
                )


if __name__ == "__main__":
    unittest.main()
