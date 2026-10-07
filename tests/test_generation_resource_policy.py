"""Unindexed fixtures must never turn into production whole-file scans."""

import gzip
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_population import indexed_lines
from generation_bindings import (
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)
from generation_resource_policy import MAX_UNINDEXED_FIXTURE_BYTES, check_unindexed_fixture
from generation_scope import prepare_generation_scope


class GenerationResourcePolicyTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / "fixture.tsv"

    def test_small_plain_fixture_can_be_read(self):
        self.path.write_text("1\t1\tA\tC\t0\t1\n")
        self.assertEqual(list(indexed_lines(self.path, None, "1", 0, 10)), [self.path.read_text()])

    def test_oversized_fallback_is_rejected_before_opening_text(self):
        with self.path.open("wb") as handle:
            handle.truncate(MAX_UNINDEXED_FIXTURE_BYTES + 1)
        with patch("generation_population.open_text") as open_text:
            with self.assertRaisesRegex(ValueError, "requires an explicit.*index"):
                list(indexed_lines(self.path, None, "1", 0, 10))
            open_text.assert_not_called()

    def test_compressed_fallback_is_rejected_even_with_plain_suffix(self):
        self.path.write_bytes(gzip.compress(b"1\t1\tA\tC\t0\t1\n"))
        with self.assertRaisesRegex(ValueError, "Compressed generation resource"):
            check_unindexed_fixture(self.path)

    def test_compressed_suffix_is_rejected_without_reading_file(self):
        for suffix in (".gz", ".bgz", ".bgzf", ".bcf"):
            with self.subTest(suffix=suffix), self.assertRaisesRegex(ValueError, "requires an explicit"):
                check_unindexed_fixture(self.path.with_suffix(suffix))

    def test_scope_rejects_oversized_fixture_before_reading_units(self):
        with self.path.open("wb") as handle:
            handle.truncate(MAX_UNINDEXED_FIXTURE_BYTES + 1)
        inputs = GenerationScopeInputs(
            population_vcf=self.path,
            population_vcf_index=None,
            score_snv=self.path,
            score_snv_index=None,
            exclusion=None,
            testing_units=self.path,
            par_regions=self.path,
            exclusions={},
            fasta=self.path,
            fasta_fai=self.path,
            mutation_rate=None,
            mutation_rate_index=None,
        )
        outputs = GenerationScopeOutputs(
            partitions=self.path,
            buckets=self.path,
            units=self.path,
            unit_intervals=self.path,
            masked_unit_intervals=self.path,
            mask=self.path,
        )
        params = GenerationScopeParameters(0.5, 10, 10, "zstd")
        with patch("generation_scope.read_bed_units") as read_units:
            with self.assertRaisesRegex(ValueError, "up to 1 MiB"):
                prepare_generation_scope(inputs, outputs, params)
            read_units.assert_not_called()


if __name__ == "__main__":
    unittest.main()
