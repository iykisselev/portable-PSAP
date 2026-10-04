"""Cross-language contracts for genomic identifiers used in bundle artifacts."""

from __future__ import annotations

import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from genomic_contracts import allele_key, allele_key_text, canonical_chromosome  # noqa: E402
from validation import canonical_chromosome_name  # noqa: E402


class GenomicContractTests(unittest.TestCase):
    def test_python_canonical_aliases_and_allele_keys(self) -> None:
        aliases = [" chr1 ", "CHR1", "x", "chrX", "M", "chrMT"]
        self.assertEqual(
            [canonical_chromosome(value) for value in aliases],
            ["1", "1", "X", "X", "MT", "MT"],
        )
        self.assertEqual(allele_key("chrM", "7", "a", "t"), ("MT", 7, "A", "T"))
        self.assertEqual(allele_key_text(("MT", 7, "A", "T")), "MT:7:A:T")

    def test_validation_applies_scope_after_shared_canonicalization(self) -> None:
        self.assertEqual(canonical_chromosome_name(" chr1 "), "1")
        self.assertEqual(canonical_chromosome_name("chrX"), "X")
        self.assertIsNone(canonical_chromosome_name("chrM"))

    def test_r_and_python_contracts_match(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")
        source = ROOT / "workflow" / "lib" / "r" / "genomic_contracts.R"
        expression = (
            f"source({source.as_posix()!r}); "
            "cat(paste(canonical_chromosome(c(' chr1 ', 'CHR1', 'x', 'chrX', 'M', 'chrMT')), collapse='\\n')); "
            "cat('\\n', canonical_allele_key('chrM', 7, 'a', 't'), sep='')"
        )
        result = subprocess.run(
            [rscript, "-e", expression], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["1", "1", "X", "X", "MT", "MT", "MT:7:A:T"])


if __name__ == "__main__":
    unittest.main()
