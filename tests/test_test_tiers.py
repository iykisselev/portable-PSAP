"""Prevent tests from silently disappearing from the CI tier selection."""

from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

from run_tier import LINT_COMMANDS, ROOT, TIERS, main, run_lint


class TestTierCoverageTests(unittest.TestCase):
    def test_every_python_test_is_assigned_once(self):
        discovered = {path.name for path in Path(__file__).parent.glob("test_*.py")}
        assigned = [name for patterns in TIERS.values() for name in patterns]
        self.assertEqual(discovered, set(assigned))
        self.assertEqual(len(assigned), len(set(assigned)))


class LocalLintTests(unittest.TestCase):
    def test_lint_runs_every_check_and_reports_failure(self):
        with patch("run_tier.subprocess.run") as run:
            run.side_effect = [
                subprocess.CompletedProcess(command, int(index == 1))
                for index, command in enumerate(LINT_COMMANDS)
            ]
            self.assertEqual(1, run_lint())
        self.assertEqual(len(LINT_COMMANDS), run.call_count)
        for command in LINT_COMMANDS:
            run.assert_any_call(command, cwd=ROOT, check=False)

    def test_missing_tool_fails_with_setup_instructions(self):
        with (
            patch("run_tier.subprocess.run", side_effect=FileNotFoundError),
            patch("sys.stderr") as stderr,
        ):
            self.assertEqual(1, run_lint())
        self.assertIn("workflow/envs/ci.yaml", str(stderr.write.call_args_list))

    def test_all_requires_lint_to_pass_before_test_discovery(self):
        for lint_result in (0, 1):
            with (
                self.subTest(lint_result=lint_result),
                patch("sys.argv", ["run_tier.py", "all"]),
                patch("run_tier.run_lint", return_value=lint_result),
                patch("run_tier.run_pattern", return_value=0) as run_tests,
            ):
                self.assertEqual(lint_result, main())
                if lint_result:
                    run_tests.assert_not_called()
                else:
                    run_tests.assert_called_once_with("test_*.py")
