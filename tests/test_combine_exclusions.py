import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


SCRIPT_DIR = Path(__file__).parents[1] / "workflow" / "lib" / "python"
sys.path.insert(0, str(SCRIPT_DIR))

from exclusions import (  # noqa: E402
    interval_overlaps_variant,
    merged_interval_index,
    vcf_record_span,
    write_projected_vcf,
)


class ExclusionOverlapTests(unittest.TestCase):
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
