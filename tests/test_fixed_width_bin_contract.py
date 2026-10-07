"""Cross-language conformance tests for the published fixed-width grid."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_contract import fixed_width_bin_count  # noqa: E402


CASES_PATH = ROOT / "tests" / "fixed_width_bin_count_cases.csv"


def cases() -> list[dict[str, str]]:
    with CASES_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class FixedWidthBinContractTests(unittest.TestCase):
    def test_python_matches_the_shared_contract_table(self) -> None:
        for case in cases():
            with self.subTest(case=case["case"]):
                self.assertEqual(
                    int(case["expected_bin_count"]),
                    fixed_width_bin_count(case["score_max"], case["bin_width"]),
                )

    def test_r_matches_the_same_contract_table(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")
        expression = """
        abort <- function(...) stop(sprintf(...), call. = FALSE)
        source(commandArgs(TRUE)[1])
        cases <- read.csv(commandArgs(TRUE)[2], stringsAsFactors = FALSE)
        values <- mapply(fixed_width_bin_count, cases$score_max, cases$bin_width)
        cat(paste(values, collapse = \"\\n\"))
        """
        result = subprocess.run(
            [
                rscript,
                "-e",
                expression,
                str(ROOT / "workflow" / "lib" / "r" / "lookup_contracts.R"),
                str(CASES_PATH),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(
            0,
            result.returncode,
            msg=f"stdout:\\n{result.stdout}\\nstderr:\\n{result.stderr}",
        )
        actual = [int(value) for value in result.stdout.splitlines()]
        expected = [int(case["expected_bin_count"]) for case in cases()]
        self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
