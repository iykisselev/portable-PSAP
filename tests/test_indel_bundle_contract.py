"""Version-three generated-bundle metadata contracts for indel nulls."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import jsonschema
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bundle_contract import is_supported_contract  # noqa: E402
from indel_bundle_contract import IndelBundleContractError, validate_indel_bundle_metadata  # noqa: E402


def component(name: str) -> dict[str, object]:
    return {"path": f"{name}.tsv.gz", "sha256": "0" * 64, "bytes": 1}


def metadata() -> dict[str, object]:
    resources = {
        name: {"identity": name, "path": f"/{name}", "sha256": "1" * 64}
        for name in ("fasta", "exclusion", "testing_units", "population_vcf", "score_snv", "score_indel")
    }
    return {
        "schema_version": 3,
        "contract": "generated_bundle_v3",
        "assembly": "GRCh38",
        "contig_convention": "canonical_without_chr_prefix",
        "testing_unit_type": "gene",
        "af_info_tags": {"overall": "AF"},
        "binning": {"algorithm": "fixed_width_decimal_v1", "score_max": 70, "bin_width": 0.05, "bin_count": 1400, "maximum_in_final_bin": True},
        "model_assumptions": {
            "roles": ["het", "hom"], "variant_classes": ["snv", "indel"],
            "unobserved_policy": "zero", "mutation_conversion": None,
            "allele_events_independent": True, "par_policy": "excluded_before_event_assignment",
            "hemizygote_role": None,
            "x_chromosome": {"af_policy": "overall_proxy", "mutation_rate_policy": "reject", "hemizygote_scope": "non_PAR_chrX_only"},
        },
        "counts": {"eligible_units": 1, "excluded_units": 0, "excluded_alleles": 0, "calibration_intervals": 1},
        "resource_identities": resources,
        "components": {
            "lookup_het": component("lookup_het"), "lookup_hom": component("lookup_hom"),
            "eligible_testing_units": component("eligible"), "calibration_regions": component("regions"),
            "excluded_alleles": component("excluded_alleles"), "excluded_units": component("excluded_units"),
            "candidate_audit": component("candidate_audit"),
        },
        "validation_artifacts": {"report": "validation.json", "marker": "validation.ok"},
        "indel_universe": {
            "candidate_universe": "score_resource", "normalization": "vcf_minimal_left_aligned_v1",
            "reference_span": "zero_based_half_open_ref_span_v1", "membership": "reference_span_overlap",
            "mask_policy": "reference_span_contained", "max_reference_length": 50, "max_alternate_length": 50,
        },
        "variant_class_counts": {
            "snv": {"retained": 10, "rejected": 1, "observed": 7, "absent": 3},
            "indel": {"retained": 2, "rejected": 1, "observed": 1, "absent": 1},
        },
    }


class IndelBundleContractTests(unittest.TestCase):
    def test_v3_schema_and_semantic_contract_accept_complete_indel_metadata(self) -> None:
        document = metadata()
        schema = yaml.safe_load((ROOT / "workflow/schemas/null_distribution_bundle.schema.yaml").read_text())
        jsonschema.Draft202012Validator(schema).validate(document)
        validate_indel_bundle_metadata(document)
        self.assertTrue(is_supported_contract("generated_bundle_v3"))

    def test_v3_mutation_fallback_requires_an_indel_rate_identity(self) -> None:
        document = metadata()
        document["model_assumptions"]["unobserved_policy"] = "mutation_rate"  # type: ignore[index]
        with self.assertRaisesRegex(IndelBundleContractError, "indel mutation-rate identity"):
            validate_indel_bundle_metadata(document)

    def test_v3_refuses_a_partial_indel_universe(self) -> None:
        document = deepcopy(metadata())
        document["indel_universe"]["membership"] = "anchor_only"  # type: ignore[index]
        with self.assertRaisesRegex(IndelBundleContractError, "unsupported indel membership"):
            validate_indel_bundle_metadata(document)


if __name__ == "__main__":
    unittest.main()
