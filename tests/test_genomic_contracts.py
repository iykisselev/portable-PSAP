"""Cross-language contracts for genomic identifiers used in bundle artifacts."""

from __future__ import annotations

import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from genomic_contracts import (  # noqa: E402
    AlleleContractError,
    allele_key,
    allele_key_text,
    canonical_chromosome,
    is_left_aligned_allele,
    is_minimal_allele,
    reference_span,
    validated_allele_key,
    variant_class,
)
from validation import canonical_chromosome_name  # noqa: E402
from bundle_contract import model_roles  # noqa: E402


class GenomicContractTests(unittest.TestCase):
    def test_r_and_python_select_the_same_enabled_roles(self):
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")
        expression = 'source(commandArgs(TRUE)[1]); cat(paste(model_roles(FALSE), collapse=","), "\\n", paste(model_roles(TRUE), collapse=","), sep="")'
        result = subprocess.run(
            [
                rscript,
                "-e",
                expression,
                str(ROOT / "workflow/lib/r/bundle_contracts.R"),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(
            [",".join(model_roles(False)), ",".join(model_roles(True))],
            result.stdout.splitlines(),
        )

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

    def test_variant_classes_and_reference_spans_follow_vcf_semantics(self) -> None:
        self.assertEqual("snv", variant_class("A", "T"))
        self.assertEqual("indel", variant_class("A", "AT"))
        self.assertEqual("indel", variant_class("AT", "A"))
        self.assertEqual("unsupported", variant_class("AT", "GC"))
        self.assertEqual("unsupported", variant_class("A", "<DEL>"))
        self.assertEqual((0, 1), reference_span(1, "A"))
        self.assertEqual((6, 8), reference_span(7, "AC"))
        self.assertTrue(is_minimal_allele("AC", "A"))
        self.assertFalse(is_minimal_allele("ATC", "AC"))
        self.assertTrue(is_left_aligned_allele(1, "A", "AA"))
        self.assertFalse(is_left_aligned_allele(2, "A", "AA"))

    def test_validated_allele_key_rejects_unsupported_alleles(self) -> None:
        self.assertEqual(
            ("1", 7, "AC", "A"),
            validated_allele_key("chr1", 7, "ac", "a"),
        )
        with self.assertRaisesRegex(AlleleContractError, "SNV or length-changing"):
            validated_allele_key("1", 7, "AC", "GT")
        with self.assertRaisesRegex(AlleleContractError, "position must be positive"):
            validated_allele_key("1", 0, "A", "T")

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
            [rscript, "-e", expression],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertEqual(
            result.stdout.splitlines(), ["1", "1", "X", "X", "MT", "MT", "MT:7:A:T"]
        )


if __name__ == "__main__":
    unittest.main()
