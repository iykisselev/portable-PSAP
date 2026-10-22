"""Regression coverage for the standalone cohort-summary review utility."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


try:
    import pysam
except ImportError:  # pragma: no cover - depends on the selected test environment
    pysam = None


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "annotate_cohort_summary.py"


@unittest.skipUnless(pysam, "pysam is required for indexed annotation tests")
class CohortSummaryAnnotationTests(unittest.TestCase):
    def test_ref_and_alt_af_cardinality_and_multiallelic_scalar_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cadd = self.tabix_table(root, "cadd.tsv", "1\t10\tA\tG\t1\t25\n")
            module = self.load_module()
            for number, frequencies, expected in (
                ("R", "0.7,0.1,0.2", "0.2"),
                ("1", "0.1", None),
            ):
                with self.subTest(number=number):
                    vcf = root / "population.vcf"
                    vcf.write_text(
                        "##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n"
                        f"##INFO=<ID=AF,Number={number},Type=Float,Description=frequency>\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        f"1\t10\t.\tA\tG,T\t.\tPASS\tAF={frequencies}\n"
                    )
                    bgzf = root / "population.vcf.gz"
                    pysam.tabix_compress(str(vcf), str(bgzf), force=True)
                    pysam.tabix_index(str(bgzf), preset="vcf", force=True)
                    lookup = module.AnnotationLookup(cadd, None, bgzf, "AF")
                    try:
                        if expected is None:
                            with self.assertRaisesRegex(ValueError, "cardinality"):
                                lookup.values(("1", 10, "A", "T"))
                        else:
                            self.assertEqual(
                                expected, lookup.values(("1", 10, "A", "T"))[1]
                            )
                    finally:
                        lookup.close()

    def test_direct_annotation_rejects_same_path_and_hardlink_without_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "summary.tsv"
            source.write_text("original content\n")
            hardlink = Path(tmp) / "hardlink.tsv"
            os.link(source, hardlink)
            module = self.load_module()
            for output in (source, hardlink):
                with self.assertRaisesRegex(ValueError, "overwrite"):
                    module.annotate_summary(
                        source, output, {}, None, "testing_unit_id", "variant_key"
                    )
                self.assertEqual("original content\n", source.read_text())

    @staticmethod
    def tabix_table(root: Path, name: str, text: str) -> Path:
        source = root / name
        source.write_text(text, encoding="utf-8")
        compressed = source.with_suffix(source.suffix + ".gz")
        pysam.tabix_compress(str(source), str(compressed), force=True)
        pysam.tabix_index(
            str(compressed), seq_col=0, start_col=1, end_col=1, force=True
        )
        return compressed

    @staticmethod
    def load_module():
        specification = importlib.util.spec_from_file_location(
            "annotate_cohort_summary", SCRIPT
        )
        assert specification and specification.loader
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module

    def test_adds_exact_cadd_af_and_gene_name_in_signature_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            cadd_snv = self.tabix_table(
                root, "cadd_snv.tsv", "1\t10\tA\tG\t1.0\t25.0\n"
            )
            cadd_indel = self.tabix_table(
                root, "cadd_indel.tsv", "1\t20\tAA\tA\t1.0\t18.0\n"
            )
            vcf = root / "population.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=1,length=100>\n"
                "##INFO=<ID=AF,Number=A,Type=Float,Description=frequency>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t10\t.\tA\tG,T\t.\tPASS\tAF=0.001,0.002\n"
                "1\t20\t.\tAA\tA\t.\tPASS\tAF=0.0003\n",
                encoding="utf-8",
            )
            population_vcf = root / "population.vcf.gz"
            pysam.tabix_compress(str(vcf), str(population_vcf), force=True)
            pysam.tabix_index(str(population_vcf), preset="vcf", force=True)
            gtf = root / "genes.gtf"
            gtf.write_text(
                "1\tsource\tgene\t1\t100\t.\t+\t.\t"
                'gene_id "ENSG000001"; gene_name "GENE_A";\n',
                encoding="utf-8",
            )
            summary = root / "cohort_candidates.tsv.gz"
            with gzip.open(summary, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "testing_unit_id\tvariant_key\tPSAP\n"
                    "ENSG000001.7\t1:10-10_A>G\t0.01\n"
                )

            module = self.load_module()
            output = root / "review.tsv.gz"
            module.main(
                [
                    str(summary),
                    str(output),
                    "--cadd-snv",
                    str(cadd_snv),
                    "--cadd-indel",
                    str(cadd_indel),
                    "--population-vcf",
                    str(population_vcf),
                    "--gtf",
                    str(gtf),
                ]
            )

            with gzip.open(output, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(1, len(rows))
            self.assertEqual("GENE_A", rows[0]["gene_name"])
            self.assertEqual("25", rows[0]["cadd_phred"])
            self.assertEqual("0.001", rows[0]["population_af"])


if __name__ == "__main__":
    unittest.main()
