"""Tests for shell fragments constructed by workflow helpers."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from shell_arguments import optional_regions_argument  # noqa: E402


class OptionalRegionsArgumentTests(unittest.TestCase):
    def test_absent_path_has_no_argument(self) -> None:
        self.assertEqual("", optional_regions_argument(None))

    def test_path_is_quoted_as_one_bcftools_argument(self) -> None:
        self.assertEqual(
            "--regions-file 'resources/calibration regions.bed'",
            optional_regions_argument(Path("resources/calibration regions.bed")),
        )


if __name__ == "__main__":
    unittest.main()
