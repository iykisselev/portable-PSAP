import gzip
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


SCRIPT_DIR = Path(__file__).parents[1] / "workflow" / "lib" / "python"
sys.path.insert(0, str(SCRIPT_DIR))

from annotation_qc import summarize_annotations  # noqa: E402


class AnnotationSummaryTests(unittest.TestCase):
    @staticmethod
    def write_annotations(path: Path, rows: list[str]) -> None:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("CHROM\tPOS\tREF\tALT\tCADD_PHRED\tunit\n")
            handle.writelines(f"{row}\n" for row in rows)

    def test_complete_snv_annotation_passes(self) -> None:
        with TemporaryDirectory() as temp_dir:
            annotations = Path(temp_dir) / "annotations.tsv.gz"
            self.write_annotations(
                annotations,
                ["chr1\t10\tA\tG\t12.4\tGENE1"],
            )

            report = summarize_annotations(annotations, include_indels=False)

            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["errors"], [])
            self.assertEqual(report["summary"]["fully_annotated_rows"], 1)

    def test_total_cadd_and_unit_miss_fails(self) -> None:
        with TemporaryDirectory() as temp_dir:
            annotations = Path(temp_dir) / "annotations.tsv.gz"
            self.write_annotations(
                annotations,
                ["chr1\t10\tA\tG\t.\t."],
            )

            report = summarize_annotations(annotations, include_indels=False)

            self.assertEqual(report["status"], "failed")
            self.assertEqual(len(report["errors"]), 2)
            self.assertTrue(any("CADD SNV" in error for error in report["errors"]))
            self.assertTrue(any("Testing-unit" in error for error in report["errors"]))

    def test_enabled_indel_annotation_must_match_an_indel(self) -> None:
        with TemporaryDirectory() as temp_dir:
            annotations = Path(temp_dir) / "annotations.tsv.gz"
            self.write_annotations(
                annotations,
                ["chr1\t10\tA\tAT\t.\tGENE1"],
            )

            report = summarize_annotations(annotations, include_indels=True)

            self.assertEqual(report["status"], "failed")
            self.assertTrue(any("CADD indel" in error for error in report["errors"]))

    def test_empty_variant_table_is_valid(self) -> None:
        with TemporaryDirectory() as temp_dir:
            annotations = Path(temp_dir) / "annotations.tsv.gz"
            self.write_annotations(annotations, [])

            report = summarize_annotations(annotations, include_indels=True)

            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["summary"]["variant_rows"], 0)


if __name__ == "__main__":
    unittest.main()
