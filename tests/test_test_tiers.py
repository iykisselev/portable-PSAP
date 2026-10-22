"""Prevent tests from silently disappearing from the CI tier selection."""

from pathlib import Path
import unittest

from run_tier import TIERS


class TestTierCoverageTests(unittest.TestCase):
    def test_every_python_test_is_assigned_once(self):
        discovered = {path.name for path in Path(__file__).parent.glob("test_*.py")}
        assigned = [name for patterns in TIERS.values() for name in patterns]
        # Partition-contract tests are also invoked directly by the focused
        # generation/apply verification command.  Keep that small suite out of
        # the historical tier manifest while still preventing any other test
        # from silently falling through the selection.
        standalone_partition_tests = {
            "test_partition_audit.py",
            "test_partition_scope.py",
            "test_partitioned_bundle_contract.py",
            "test_partitioned_lookup.py",
        }
        self.assertEqual(discovered - standalone_partition_tests, set(assigned))
        self.assertEqual(len(assigned), len(set(assigned)))
