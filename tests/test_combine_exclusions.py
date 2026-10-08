import sys
import unittest
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch


SCRIPT_DIR = Path(__file__).parents[1] / "workflow" / "lib" / "python"
sys.path.insert(0, str(SCRIPT_DIR))

import exclusions  # noqa: E402
from exclusions import (  # noqa: E402
    collect_exclusions,
    interval_overlaps_variant,
    merged_interval_index,
    project_vcf_records,
    vcf_record_span,
    write_projected_vcf,
)
from validation import validate_exclusion_tsv  # noqa: E402


class ExclusionOverlapTests(unittest.TestCase):
    def test_validation_and_staging_share_tsv_target_semantics(self) -> None:
        with TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "targets.tsv"
            source.write_text(
                "# comment\ntrack name=ignored\nbrowser position ignored\n"
                "chr1\t10\t12\nchrM\t15\tA\tC,T\n",
                encoding="utf-8",
            )
            errors: list[str] = []
            chromosomes = validate_exclusion_tsv(source, "Exclusion TSV", errors)
            intervals, variants, indexed_sources = collect_exclusions([f"tsv:{source}"])

        self.assertEqual([], errors)
        self.assertEqual({"1", "MT"}, chromosomes)
        self.assertEqual(([10], [12]), intervals["1"])
        self.assertEqual(
            {("MT", 15, "A", "C"), ("MT", 15, "A", "T")}, variants
        )
        self.assertEqual((), indexed_sources)

    def test_validation_and_staging_reject_the_same_bad_tsv_row(self) -> None:
        with TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "targets.tsv"
            source.write_text("1\t4\tA\tC,\n", encoding="utf-8")
            errors: list[str] = []
            validate_exclusion_tsv(source, "Exclusion TSV", errors)
            with self.assertRaisesRegex(ValueError, "empty ALT allele"):
                collect_exclusions([f"tsv:{source}"])

        self.assertEqual(["Exclusion TSV line 1 has an empty ALT allele"], errors)

    def test_large_indexed_tsv_is_queried_without_building_an_allele_set(self) -> None:
        with TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "excluded_alleles.tsv.gz"
            Path(f"{source}.tbi").touch()
            intervals, variants, indexed_sources = collect_exclusions(
                [f"tsv:{source}"], indexed_tsv_paths=[str(source)]
            )
            reader = MagicMock()
            reader.fetch.return_value = iter(["1\t10\tA\tC\tfiltered\n"])
            tabix = MagicMock()
            tabix.return_value.__enter__.return_value = reader
            output = StringIO()
            with patch.object(exclusions, "pysam", MagicMock(TabixFile=tabix)):
                project_vcf_records(
                    [
                        "##fileformat=VCFv4.3\n",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                        "1\t10\t.\tA\tC\t.\tPASS\t.\n",
                        "1\t11\t.\tA\tG\t.\tPASS\t.\n",
                    ],
                    output,
                    intervals,
                    variants,
                    "cohort.vcf",
                    indexed_tsv_paths=indexed_sources,
                )

        self.assertFalse(variants)
        self.assertEqual((source,), indexed_sources)
        self.assertEqual(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t10\t.\tA\tC\t.\tPASS\t.\n",
            output.getvalue(),
        )
        reader.fetch.assert_any_call("1", 9, 10)

    def test_deletion_span_overlaps_interval_after_its_position(self) -> None:
        index = merged_interval_index({"1": [(100, 110)]})

        self.assertTrue(interval_overlaps_variant(index, "chr1", 95, 104))
        self.assertFalse(interval_overlaps_variant(index, "chr1", 95, 99))
        self.assertFalse(interval_overlaps_variant(index, "chr1", 111, 120))

    def test_mitochondrial_aliases_share_the_genomic_contract(self) -> None:
        index = merged_interval_index({"MT": [(1, 10)]})
        self.assertTrue(interval_overlaps_variant(index, "chrM", 5, 5))

    def test_vcf_span_uses_ref_length_and_info_end(self) -> None:
        self.assertEqual(
            vcf_record_span(["1", "95", ".", "AAAAAAAAAA", "A"]),
            (95, 104),
        )
        self.assertEqual(
            vcf_record_span(
                ["1", "90", ".", "N", "<DEL>", ".", "PASS", "END=105"]
            ),
            (90, 105),
        )

    def test_projected_vcf_includes_full_span_overlap(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cohort = temp_path / "cohort.vcf"
            output = temp_path / "excluded.vcf"
            cohort.write_text(
                "##fileformat=VCFv4.3\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t95\t.\tAAAAAAAAAA\tA\t.\tPASS\t.\n"
                "1\t111\t.\tA\tG\t.\tPASS\t.\n",
                encoding="utf-8",
            )

            write_projected_vcf(
                cohort,
                output,
                merged_interval_index({"1": [(100, 110)]}),
                set(),
            )

            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "##fileformat=VCFv4.3\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t95\t.\tAAAAAAAAAA\tA\t.\tPASS\t.\n",
            )


if __name__ == "__main__":
    unittest.main()
