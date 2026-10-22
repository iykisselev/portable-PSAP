import gzip
import shutil
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch


SCRIPT_DIR = Path(__file__).parents[1] / "workflow" / "lib" / "python"
sys.path.insert(0, str(SCRIPT_DIR))

from validation import (  # noqa: E402
    IndexedResource,
    VcfHeaderFacts,
    read_vcf_header_facts,
    read_ped,
    read_samples_tsv,
    validate_cadd_chromosome_names,
    validate_fasta,
    validate_hemizygote_chromosomes,
    validate_indexed_table,
    validate_reference_tables,
    validate_region_chromosomes,
)
import resources  # noqa: E402
from validation_cohort import validate_cohort_inputs  # noqa: E402
from validation_resources import validate_resource_inputs  # noqa: E402


class ValidationTests(unittest.TestCase):
    def test_vcf_facts_use_index_statistics_when_contig_headers_are_absent(
        self,
    ) -> None:
        with TemporaryDirectory() as temp_dir:
            vcf = Path(temp_dir) / "cohort.vcf.gz"
            vcf.write_bytes(b"placeholder")
            header = SimpleNamespace(
                returncode=0,
                stdout="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
                stderr="",
            )
            statistics = SimpleNamespace(
                returncode=0,
                stdout="1\t100\nX\t50\n",
                stderr="",
            )
            errors: list[str] = []
            with (
                patch("validation.shutil.which", return_value="bcftools"),
                patch(
                    "validation.subprocess.run", side_effect=(header, statistics)
                ) as run,
                patch("validation.read_vcf_record_contigs") as record_contigs,
            ):
                facts = read_vcf_header_facts(vcf, errors)

        self.assertEqual(errors, [])
        self.assertEqual(facts.samples, ("sample",))
        self.assertEqual(facts.contigs, frozenset({"1", "X"}))
        self.assertEqual(facts.canonical_map, {"1": "1", "X": "X"})
        self.assertEqual(
            run.call_args_list[1].args[0], ["bcftools", "index", "--stats", str(vcf)]
        )
        record_contigs.assert_not_called()

    def test_all_unknown_sample_sheet_phenotypes_fail_validation(self) -> None:
        with TemporaryDirectory() as temp_dir:
            samples = Path(temp_dir) / "samples.tsv"
            samples.write_text(
                "family_id\tpid\tsex\tphenotype\n"
                "F1\tunknown1\t1\t0\n"
                "F1\tunknown2\t2\t-9\n",
                encoding="utf-8",
            )
            errors: list[str] = []
            warnings: list[str] = []

            result = read_samples_tsv(samples, errors, warnings)

            self.assertEqual(result["sample_ids"], [])
            self.assertEqual(warnings, [])
            self.assertEqual(
                errors,
                [
                    "Sample sheet contains no samples with a known phenotype "
                    "(1=unaffected or 2=affected)"
                ],
            )

    def test_all_unknown_phenotypes_fail_validation(self) -> None:
        with TemporaryDirectory() as temp_dir:
            ped = Path(temp_dir) / "cohort.ped"
            ped.write_text(
                "F1\tunknown1\t0\t0\t1\t0\nF1\tunknown2\t0\t0\t2\t-9\n",
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

    def test_hemizygote_model_requires_x(self) -> None:
        errors: list[str] = []

        validate_hemizygote_chromosomes(True, {"1": "1"}, errors)

        self.assertEqual(
            errors,
            ["HEM chromosome-X model requires chromosome(s) in the input VCF: X"],
        )

    def test_sex_chromosomes_are_optional_when_hemizygote_model_is_disabled(
        self,
    ) -> None:
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

    def test_fasta_uses_the_declared_index_binding(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "reference.fa"
            declared_index = temp_path / "reference.index.fai"
            fasta.write_text(">chr1\nA\n", encoding="utf-8")
            declared_index.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")
            errors: list[str] = []

            validate_fasta(fasta, declared_index, None, errors)

            self.assertEqual(errors, [])
            self.assertFalse((temp_path / "reference.fa.fai").exists())

    def test_table_uses_the_declared_index_binding(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            table = temp_path / "resource.tsv"
            declared_index = temp_path / "resource.index.tbi"
            table.write_text("chr1\t10\t20\n", encoding="utf-8")
            declared_index.write_text("index", encoding="utf-8")
            errors: list[str] = []
            warnings: list[str] = []

            validate_indexed_table(
                table,
                "Test table",
                errors,
                warnings,
                declared_index,
            )

            self.assertEqual(errors, [])
            self.assertFalse(
                any("no configured index" in warning for warning in warnings)
            )
            self.assertFalse((temp_path / "resource.tsv.tbi").exists())

    def test_reference_tables_preserve_comma_paths_and_index_bindings(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            unindexed_exclusion = temp_path / "first,exclusion.bed"
            indexed_exclusion = temp_path / "second-exclusion.bed"
            exclusion_index = temp_path / "second-exclusion.index.tbi"
            par_regions = temp_path / "par.bed"
            for path, contents in (
                (unindexed_exclusion, "1\t0\t5\n"),
                (indexed_exclusion, "1\t5\t10\n"),
                (par_regions, "X\t0\t10\n"),
            ):
                path.write_text(contents, encoding="utf-8")
            exclusion_index.write_text("index\n", encoding="utf-8")
            errors: list[str] = []
            warnings: list[str] = []

            chromosomes = validate_reference_tables(
                calibration_regions=None,
                exclusions=(
                    IndexedResource(unindexed_exclusion),
                    IndexedResource(indexed_exclusion, exclusion_index),
                ),
                par_regions=(IndexedResource(par_regions),),
                errors=errors,
                warnings=warnings,
            )

            self.assertEqual(errors, [])
            self.assertEqual(chromosomes["exclusions"], {"1", "X"})
            self.assertTrue(
                any("first,exclusion.bed" in warning for warning in warnings)
            )
            self.assertTrue(
                any(
                    "Exclusion file 1 has no configured index" in warning
                    for warning in warnings
                )
            )
            self.assertFalse(
                any(
                    "Exclusion file 2 has no configured index" in warning
                    for warning in warnings
                )
            )
            self.assertTrue(
                any(
                    "Exclusion file 3 has no configured index" in warning
                    for warning in warnings
                )
            )

    def test_resource_validation_forwards_indexed_resource_bindings(self) -> None:
        exclusions = (
            IndexedResource("first,exclusion.bed"),
            IndexedResource("second-exclusion.bed", "second-exclusion.tbi"),
        )
        with (
            patch(
                "validation_resources.validate_reference_tables",
                return_value={"exclusions": {"1"}},
            ) as validate_tables,
            patch(
                "validation_resources.read_vcf_header_facts",
                return_value=VcfHeaderFacts((), frozenset({"1"}), {"1": "1"}),
            ),
            patch(
                "validation_resources.validate_cadd_resources", return_value={}
            ) as validate_cadd,
        ):
            report = validate_resource_inputs(
                vcf="cohort.vcf.gz",
                calibration_regions=None,
                exclusions=exclusions,
                par_regions=(),
                cadd_snv=IndexedResource("cadd.tsv.gz", "cadd.tsv.gz.tbi"),
                cadd_indel=None,
                include_indels=False,
            )

        self.assertEqual(validate_tables.call_args.kwargs["exclusions"], exclusions)
        self.assertIs(validate_cadd.call_args.args[2], False)
        self.assertIs(report["summary"]["include_indels"], False)

    def test_cohort_validation_passes_a_native_hemizygote_boolean(self) -> None:
        metadata = {
            "manifest": [],
            "all_sample_ids": [],
            "excluded_phenotype_sample_ids": [],
            "sex_counts": {},
            "phenotype_counts": {},
        }
        facts = VcfHeaderFacts((), frozenset({"1"}), {"1": "1"})
        with (
            patch("validation_cohort.read_samples_tsv", return_value=metadata),
            patch("validation_cohort.write_sample_manifest"),
            patch(
                "validation_cohort.read_sample_manifest_index",
                return_value=SimpleNamespace(records=[]),
            ),
            patch(
                "validation_cohort.read_vcf_header_facts",
                return_value=facts,
            ),
            patch(
                "validation_cohort.validate_hemizygote_chromosomes"
            ) as validate_hemizygote,
            patch("validation_cohort.validate_vcf_samples"),
            patch("validation_cohort.validate_vcf_resource"),
            patch("validation_cohort.validate_fasta"),
            patch(
                "validation_cohort.validate_normalized_split_vcf"
            ) as validate_normalized,
        ):
            validate_cohort_inputs(
                metadata_path="samples.tsv",
                metadata_kind="samples",
                vcf="cohort.vcf.gz",
                vcf_index=None,
                fasta="reference.fa",
                fasta_index="reference.fa.fai",
                fasta_gzi=None,
                hemizygote_model=True,
                manifest_path="sample_manifest.tsv",
            )

        self.assertIs(validate_hemizygote.call_args.args[0], True)
        self.assertIs(validate_normalized.call_args.args[-1], facts)

    def test_staging_accepts_a_declared_nonconventional_index(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            declared_index = temp_path / "resource.index.tbi"
            declared_index.write_text("index", encoding="utf-8")

            self.assertTrue(resources.source_has_index(declared_index))
            self.assertFalse(
                resources.is_conventional_source_index(
                    temp_path / "resource.tsv.gz", declared_index
                )
            )

    def test_staging_reuses_a_complete_indexed_vcf_without_resorting(self) -> None:
        with TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            source = temp_path / "cohort.vcf.gz"
            source_index = Path(f"{source}.tbi")
            output = temp_path / "prepared.vcf.gz"
            source.write_bytes(b"placeholder")
            source_index.write_text("index", encoding="utf-8")

            with (
                patch.object(resources, "is_bgzf", return_value=True),
                patch.object(resources, "has_bgzf_eof", return_value=True),
                patch.object(resources, "require_tool", return_value="bcftools"),
                patch.object(resources, "run") as run,
            ):
                resources.prepare_vcf(source, output, source_index)

            self.assertTrue(output.is_symlink())
            self.assertEqual(source.resolve(), output.resolve())
            self.assertEqual(
                ["bcftools", "index", "--force", "--tbi", str(output)],
                run.call_args.args[0],
            )

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
                resources,
                "require_tool",
                side_effect=lambda name: (
                    shutil.which("gzip") if name == "bgzip" else "/usr/bin/sort"
                ),
            ):
                resources.copy_as_bgzip(source, output)

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
