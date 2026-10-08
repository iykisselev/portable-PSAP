"""Cross-language conformance for missing table-cell and identifier tokens."""

from __future__ import annotations

import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from io_contracts import (  # noqa: E402
    MISSING_TOKENS,
    PED_MISSING_SEX_TOKENS,
    is_missing_token,
)


EXPECTED_MISSING_TOKENS = {"", ".", "-", "na", "nan", "none", "null"}
EXPECTED_PED_MISSING_SEX_TOKENS = EXPECTED_MISSING_TOKENS | {"-9", "0"}


def decode_r_tokens(output: str) -> tuple[set[str], set[str]]:
    """Decode length-prefixed R strings, preserving the empty token."""
    ordinary, separator, ped = output.partition("--PED--\n")
    if not separator:
        raise AssertionError(f"R token output has no PED separator: {output!r}")

    def decode(lines: str) -> set[str]:
        tokens = set()
        for line in lines.splitlines():
            length_text, value = line.split(":", 1)
            if int(length_text) != len(value):
                raise AssertionError(f"Malformed length-prefixed R token: {line!r}")
            tokens.add(value)
        return tokens

    return decode(ordinary), decode(ped)


class MissingTokenContractTests(unittest.TestCase):
    def test_python_tokens_and_ped_extension_are_deliberately_distinct(self) -> None:
        self.assertEqual(EXPECTED_MISSING_TOKENS, MISSING_TOKENS)
        self.assertEqual(EXPECTED_PED_MISSING_SEX_TOKENS, PED_MISSING_SEX_TOKENS)
        for token in EXPECTED_MISSING_TOKENS:
            self.assertTrue(is_missing_token(token.upper()))
        self.assertFalse(is_missing_token("0"))
        self.assertFalse(is_missing_token("-9"))

    def test_r_uses_the_same_token_sets(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")
        expression = """
        abort <- function(...) stop(sprintf(...), call. = FALSE)
        source(commandArgs(TRUE)[1])
        emit <- function(tokens) cat(paste(sprintf("%d:%s", nchar(tokens), tokens), collapse = "\\n"), "\\n", sep = "")
        emit(MISSING_TOKENS)
        cat("--PED--\\n")
        emit(PED_MISSING_SEX_TOKENS)
        """
        result = subprocess.run(
            [
                rscript,
                "-e",
                expression,
                str(ROOT / "workflow" / "lib" / "r" / "table_contracts.R"),
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
        r_missing, r_ped_sex = decode_r_tokens(result.stdout)
        self.assertEqual(EXPECTED_MISSING_TOKENS, r_missing)
        self.assertEqual(EXPECTED_PED_MISSING_SEX_TOKENS, r_ped_sex)


if __name__ == "__main__":
    unittest.main()
