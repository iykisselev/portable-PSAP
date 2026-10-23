"""Maintenance contracts for documented run-mode examples."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from workflow_plan import build_workflow_plan  # noqa: E402


def merge(base: dict, override: dict) -> dict:
    """Recursively apply one run manifest over the scientific defaults."""
    result = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


class StageEightContracts(unittest.TestCase):
    def test_every_supported_mode_has_a_current_example_plan(self) -> None:
        science = yaml.safe_load(
            (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
        )
        examples = {
            "generate": ROOT / "runs" / "generate.example.yaml",
            "full": ROOT / "runs" / "full.example.yaml",
        }
        for mode, path in examples.items():
            with self.subTest(mode=mode):
                plan = build_workflow_plan(
                    merge(science, yaml.safe_load(path.read_text(encoding="utf-8")))
                )
                self.assertEqual(mode, plan.mode)
        apply_example = yaml.safe_load(
            (ROOT / "runs" / "apply.example.yaml").read_text(encoding="utf-8")
        )
        self.assertEqual(
            "bundle", apply_example["bindings"]["reference"]["lookup_source"]["type"]
        )


if __name__ == "__main__":
    unittest.main()
