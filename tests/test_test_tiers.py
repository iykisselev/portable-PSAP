"""Prevent tests from silently disappearing from the CI tier selection."""

from pathlib import Path
import unittest

from run_tier import TIERS


class TestTierCoverageTests(unittest.TestCase):
    def test_every_python_test_is_assigned_once(self):
        discovered = {path.name for path in Path(__file__).parent.glob("test_*.py")}
        assigned = [name for patterns in TIERS.values() for name in patterns]
        self.assertEqual(discovered, set(assigned))
        self.assertEqual(len(assigned), len(set(assigned)))
