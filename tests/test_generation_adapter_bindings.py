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
from apply_bindings import (
    CohortValidationInputs,
    PreparationInputs,
    AnnotationParameters,
)  # noqa: E402
from generation_bindings import (  # noqa: E402
    BundleFinalizationParameters,
    GenerationScopeInputs,
    GenerationScopeParameters,
    LookupReductionInputs,
    LookupReductionOutputs,
)


class GenerationAdapterBindingTests(unittest.TestCase):
    def test_boolean_binding_rejects_null_and_arbitrary_truthy_values(self):
        for value in (None, 2, [], ["false"]):
            with (
                self.subTest(value=value),
                self.assertRaisesRegex(BindingError, "include_indels"),
            ):
                bind(AnnotationParameters, SimpleNamespace(include_indels=value))

    def test_optional_bindings_reject_extra_paths_with_field_context(self):
        with self.assertRaisesRegex(
            BindingError, "PreparationInputs.source_index.*at most one"
        ):
            bind(
                PreparationInputs,
                SimpleNamespace(source="x", source_index=["one", "two"]),
            )

    def test_invalid_scalar_binding_reports_the_field(self):
        from generation_bindings import LookupReductionParameters

        with self.assertRaisesRegex(
            BindingError, "LookupReductionParameters.bin_count"
        ):
            bind(
                LookupReductionParameters,
                SimpleNamespace(bin_count="not-integer", roles=["het"]),
            )

    def test_apply_inputs_decode_optional_paths_and_require_core_fields(self):
        preparation = bind(
            PreparationInputs, SimpleNamespace(source="resource.vcf", source_index=[])
        )
        self.assertIsNone(preparation.source_index)
        self.assertIsNone(preparation.cohort)
        cohort = bind(
            CohortValidationInputs,
            SimpleNamespace(
                metadata="cohort.ped",
                vcf="cohort.vcf",
                vcf_index=["custom.tbi"],
                fasta="reference.fa",
                fasta_index="reference.fa.fai",
                fasta_gzi=[],
            ),
        )
        self.assertEqual("custom.tbi", cohort.vcf_index)
        self.assertIsNone(cohort.fasta_gzi)
        with self.assertRaisesRegex(BindingError, "source"):
            bind(PreparationInputs, SimpleNamespace(source_index=[]))

    def test_scope_binding_normalizes_paths_optionals_mappings_and_scalars(
        self,
    ) -> None:
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
                exclusions=[Path("low-complexity.bed"), "blacklist.bed"],
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
        self.assertEqual(("low-complexity.bed", "blacklist.bed"), inputs.exclusions)
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
                lookup_hem_chrX=[],
                diagnostics="diagnostics.parquet",
            ),
        )

        self.assertEqual(("one.parquet", "two.parquet"), inputs.partition_statistics)
        self.assertIsNone(outputs.lookup_hem_chrX)

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
                resource_identities={
                    "fasta": {"identity": Path("GRCh38"), "path": Path("reference.fa")}
                },
                indel_candidate_universe="score_resource",
                indel_membership="reference_span_overlap",
                indel_mask_policy="reference_span_contained",
                indel_max_reference_length=50,
                indel_max_alternate_length=50,
            ),
        )

        self.assertEqual(("snv", "indel"), parameters.variant_classes)
        self.assertEqual({"AC": "AC"}, parameters.af_tags)
        self.assertEqual(
            {"fasta": {"identity": "GRCh38", "path": "reference.fa"}},
            parameters.resource_identities,
        )
        self.assertTrue(parameters.hemizygote)
        self.assertEqual(50, parameters.indel_max_reference_length)

    def test_missing_indel_limit_fails_at_boundary(self) -> None:
        parameters = SimpleNamespace(
            minimum_retained_fraction=0.5,
            partition_size=10,
            unit_bucket_size=10,
            compression="zstd",
            af_tags={"overall": "AF"},
            include_indels=True,
        )
        with self.assertRaisesRegex(BindingError, "indel_max_reference_length"):
            bind(GenerationScopeParameters, parameters)

    def test_missing_required_binding_is_clear(self) -> None:
        with self.assertRaisesRegex(
            BindingError,
            "GenerationScopeInputs requires Snakemake field 'population_vcf'",
        ):
            bind(GenerationScopeInputs, SimpleNamespace())


if __name__ == "__main__":
    unittest.main()
