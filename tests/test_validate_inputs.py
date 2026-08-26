import gzip
import shutil
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch


SCRIPT_DIR = Path(__file__).parents[1] / "workflow" / "scripts" / "apply"
sys.path.insert(0, str(SCRIPT_DIR))

from validate_inputs import (  # noqa: E402
    read_ped,
    validate_cadd_chromosome_names,
    validate_hemizygote_chromosomes,
    validate_lookup_unit_overlap,
    validate_region_chromosomes,
)
import prepare_genomic_resource  # noqa: E402


class ValidationTests(unittest.TestCase):
    def test_samples_with_unknown_phenotype_are_excluded(self) -> None:
        with TemporaryDirectory() as temp_dir:
            ped = Path(temp_dir) / "cohort.ped"
            ped.write_text(
                "F1\taffected\t0\t0\t1\t2\n"
                "F1\tunaffected\t0\t0\t2\t1\n"
                "F1\tunknown\t0\t0\t1\t0\n",
                encoding="utf-8",
            )
            errors: list[str] = []
            warnings: list[str] = []

            result = read_ped(ped, errors, warnings)

            self.assertEqual(errors, [])
            self.assertEqual(
                result["all_sample_ids"],
                ["affected", "unaffected", "unknown"],
            )
            self.assertEqual(result["sample_ids"], ["affected", "unaffected"])
            self.assertEqual(
                [row["pid"] for row in result["manifest"]],
                ["affected", "unaffected"],
            )
            self.assertEqual(result["excluded_phenotype_sample_ids"], ["unknown"])
            self.assertEqual(len(warnings), 1)
            self.assertIn("excluded", warnings[0].lower())
            self.assertIn("unknown", warnings[0])

    def test_all_unknown_phenotypes_fail_validation(self) -> None:
        with TemporaryDirectory() as temp_dir:
            ped = Path(temp_dir) / "cohort.ped"
            ped.write_text(
                "F1\tunknown1\t0\t0\t1\t0\n"
                "F1\tunknown2\t0\t0\t2\t-9\n",
                encoding="utf-8",
            )
            errors: list[str] = []
            warnings: list[str] = []

            result = read_ped(ped, errors, warnings)

            self.assertEqual(result["sample_ids"], [])
            self.assertEqual(result["manifest"], [])
            self.assertEqual(warnings, [])
            self.assertEqual(
                errors,
                [
                    "PED file contains no samples with a known phenotype "
                    "(1=unaffected or 2=affected)"
                ],
            )

    def test_exclusion_chromosome_mismatch_is_a_warning(self) -> None:
        warnings: list[str] = []

        validate_region_chromosomes(
            {"X", "Y"},
            {"1": "1"},
            "Exclusion interval files",
            warnings,
        )

        self.assertEqual(
            warnings,
            [
                "Exclusion interval files contains chromosome name(s) absent "
                "from the input VCF: X, Y"
            ],
        )

    def test_hemizygote_model_requires_x_and_y(self) -> None:
        errors: list[str] = []

        validate_hemizygote_chromosomes(True, {"1": "1"}, errors)

        self.assertEqual(
            errors,
            [
                "Hemizygote model requires sex chromosome(s) in the input VCF: X, Y"
            ],
        )

    def test_sex_chromosomes_are_optional_when_hemizygote_model_is_disabled(self) -> None:
        errors: list[str] = []

        validate_hemizygote_chromosomes(False, {"1": "1"}, errors)

        self.assertEqual(errors, [])

    def test_cadd_alias_mismatch_is_harmonized(self) -> None:
        errors: list[str] = []
        warnings: list[str] = []

        validate_cadd_chromosome_names(
            {"1", "2"},
            {"1": "chr1", "2": "chr2"},
            "CADD SNVs",
            errors,
            warnings,
        )

        self.assertEqual(errors, [])
        self.assertEqual(len(warnings), 1)
        self.assertIn("will be harmonized", warnings[0])

    def test_cadd_without_canonical_vcf_overlap_fails(self) -> None:
        errors: list[str] = []
        warnings: list[str] = []

        validate_cadd_chromosome_names(
            {"chr1"},
            {"X": "chrX"},
            "CADD SNVs",
            errors,
            warnings,
        )

        self.assertEqual(len(errors), 1)
        self.assertIn("no canonical chromosome overlap", errors[0])
        self.assertEqual(warnings, [])

    def test_lookup_with_zero_testing_unit_overlap_fails(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bed = temp_path / "units.bed"
            lookup = temp_path / "lookup.tsv"
            bed.write_text("chr1\t0\t10\tGENE1\n", encoding="utf-8")
            lookup.write_text("OTHER\t0.1\t0.2\n", encoding="utf-8")
            errors: list[str] = []

            summary = validate_lookup_unit_overlap(
                bed,
                {"het": str(lookup), "hom": str(lookup), "hem_chrX": ""},
                {"hemizygote": False, "include_indels": False},
                errors,
            )

            self.assertEqual(summary["het"]["overlap_count"], 0)
            self.assertEqual(summary["hom"]["overlap_count"], 0)
            self.assertEqual(len(errors), 2)
            self.assertTrue(all("zero testing-unit ID overlap" in error for error in errors))

    def test_preparing_a_table_sorts_rows_before_bgzip(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            source = temp_path / "par.bed"
            output = temp_path / "par.bed.gz"
            source.write_text(
                "track name=par\n"
                "chr1\t56887904\t56887905\n"
                "chr1\t7217417\t7217418\n"
                "chr2\t1\t2\n",
                encoding="utf-8",
            )

            with patch.object(
                prepare_genomic_resource,
                "require_tool",
                side_effect=lambda name: (
                    shutil.which("gzip") if name == "bgzip" else "/usr/bin/sort"
                ),
            ):
                prepare_genomic_resource.copy_as_bgzip(source, output)

            with gzip.open(output, "rt", encoding="utf-8") as prepared:
                self.assertEqual(
                    prepared.read(),
                    "track name=par\n"
                    "chr1\t7217417\t7217418\n"
                    "chr1\t56887904\t56887905\n"
                    "chr2\t1\t2\n",
                )


if __name__ == "__main__":
    unittest.main()
