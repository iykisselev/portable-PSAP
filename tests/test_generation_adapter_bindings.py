"""Contracts for typed Snakemake adapter binding."""

from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "shared"))
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from _library import BindingError, bind  # noqa: E402
from generation_bindings import (  # noqa: E402
    BundleFinalizationParameters,
    GenerationScopeInputs,
    GenerationScopeParameters,
    LookupReductionInputs,
    LookupReductionOutputs,
)


class GenerationAdapterBindingTests(unittest.TestCase):
    def test_scope_binding_normalizes_paths_optionals_mappings_and_scalars(self) -> None:
        inputs = bind(
            GenerationScopeInputs,
            SimpleNamespace(
                population_vcf=Path("population.vcf.gz"),
                population_vcf_index=[Path("population.vcf.gz.tbi")],
                score_snv="scores.tsv.gz",
                score_snv_index=[],
                exclusion=[],
                testing_units="units.bed",
                par_regions="par.bed",
                exclusions={"low_complexity": Path("low-complexity.bed")},
                fasta="reference.fa",
                fasta_fai="reference.fa.fai",
                mutation_rate=[],
                mutation_rate_index=[],
            ),
        )
        parameters = bind(
            GenerationScopeParameters,
            SimpleNamespace(
                minimum_retained_fraction="0.75",
                partition_size="5000000",
                unit_bucket_size="100",
                compression=Path("zstd"),
                af_tags={"AC": "AC"},
                include_indels="false",
                indel_max_reference_length="12",
            ),
        )

        self.assertEqual("population.vcf.gz", inputs.population_vcf)
        self.assertEqual("population.vcf.gz.tbi", inputs.population_vcf_index)
        self.assertIsNone(inputs.score_snv_index)
        self.assertIsNone(inputs.exclusion)
        self.assertEqual({"low_complexity": "low-complexity.bed"}, inputs.exclusions)
        self.assertEqual(0.75, parameters.minimum_retained_fraction)
        self.assertEqual(5_000_000, parameters.partition_size)
        self.assertFalse(parameters.include_indels)
        self.assertEqual(12, parameters.indel_max_reference_length)

    def test_tuple_paths_and_optional_output_are_normalized(self) -> None:
        inputs = bind(
            LookupReductionInputs,
            SimpleNamespace(
                buckets=Path("buckets.json"),
                partition_statistics=[Path("one.parquet"), "two.parquet"],
            ),
        )
        outputs = bind(
            LookupReductionOutputs,
            SimpleNamespace(
                lookup_het="het.parquet",
                lookup_hom="hom.parquet",
                lookup_hem=[],
                diagnostics="diagnostics.parquet",
            ),
        )

        self.assertEqual(("one.parquet", "two.parquet"), inputs.partition_statistics)
        self.assertIsNone(outputs.lookup_hem)

    def test_mapping_and_tuple_parameters_are_typed(self) -> None:
        parameters = bind(
            BundleFinalizationParameters,
            SimpleNamespace(
                bin_count="10",
                minimum_retained_fraction="0.9",
                bin_width="0.5",
                score_max="50",
                policy="all",
                conversion="hg38",
                assembly="GRCh38",
                testing_unit_type="gene",
                af_tags={"AC": Path("AC")},
                x_af_policy="pooled",
                x_mutation_rate_policy="mean",
                hemizygote="true",
                variant_classes=[Path("snv"), "indel"],
                resource_ids={"fasta": Path("sha256")},
            ),
        )

        self.assertEqual(("snv", "indel"), parameters.variant_classes)
        self.assertEqual({"AC": "AC"}, parameters.af_tags)
        self.assertEqual({"fasta": "sha256"}, parameters.resource_ids)
        self.assertTrue(parameters.hemizygote)
        self.assertEqual(1, parameters.indel_max_reference_length)

    def test_missing_required_binding_is_clear(self) -> None:
        with self.assertRaisesRegex(
            BindingError,
            "GenerationScopeInputs requires Snakemake field 'population_vcf'",
        ):
            bind(GenerationScopeInputs, SimpleNamespace())


if __name__ == "__main__":
    unittest.main()
