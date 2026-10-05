"""Semantic requirements unique to certified combined SNV/indel bundles."""

from __future__ import annotations

from typing import Mapping


class IndelBundleContractError(ValueError):
    """Raised when generated_bundle_v3 lacks an indel scientific contract."""


def validate_indel_bundle_metadata(metadata: Mapping[str, object]) -> None:
    """Check cross-field v3 invariants beyond JSON Schema structure."""
    assumptions = metadata.get("model_assumptions")
    universe = metadata.get("indel_universe")
    counts = metadata.get("variant_class_counts")
    resources = metadata.get("resource_identities")
    components = metadata.get("components")
    if not all(isinstance(value, Mapping) for value in (assumptions, universe, counts, resources, components)):
        raise IndelBundleContractError("generated_bundle_v3 has malformed indel metadata")
    if assumptions.get("variant_classes") != ["snv", "indel"]:
        raise IndelBundleContractError("generated_bundle_v3 must certify variant_classes [snv, indel]")
    required_universe = {
        "candidate_universe": "score_resource",
        "normalization": "vcf_minimal_left_aligned_v1",
        "reference_span": "zero_based_half_open_ref_span_v1",
        "membership": "reference_span_overlap",
        "mask_policy": "reference_span_contained",
    }
    for key, expected in required_universe.items():
        if universe.get(key) != expected:
            raise IndelBundleContractError(f"generated_bundle_v3 has unsupported indel {key}")
    for key in ("max_reference_length", "max_alternate_length"):
        if not isinstance(universe.get(key), int) or universe[key] < 1:
            raise IndelBundleContractError(f"generated_bundle_v3 has invalid indel {key}")
    for variant_class in ("snv", "indel"):
        row = counts.get(variant_class)
        if not isinstance(row, Mapping) or any(
            not isinstance(row.get(name), int) or row[name] < 0
            for name in ("retained", "rejected", "observed", "absent")
        ):
            raise IndelBundleContractError(
                f"generated_bundle_v3 has invalid {variant_class} candidate counts"
            )
    if not {"score_snv", "score_indel"}.issubset(resources):
        raise IndelBundleContractError("generated_bundle_v3 lacks SNV or indel score identities")
    if assumptions.get("unobserved_policy") == "mutation_rate" and "mutation_rate_indel" not in resources:
        raise IndelBundleContractError("generated_bundle_v3 lacks the indel mutation-rate identity")
    if "candidate_audit" not in components:
        raise IndelBundleContractError("generated_bundle_v3 lacks its candidate audit component")
