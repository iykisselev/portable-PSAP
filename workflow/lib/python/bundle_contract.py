"""Versioned generated-bundle names, roles, and upgrade policy.

Only :data:`CURRENT_GENERATED_BUNDLE_CONTRACT` may be emitted by this version
of the workflow.  Older contracts are accepted by consumers through the
explicit compatibility registry; a future writer must add an upgrader rather
than silently reinterpret an existing certificate.
"""

from __future__ import annotations

import math
import sys
from types import MappingProxyType
from typing import Mapping


CURRENT_GENERATED_BUNDLE_CONTRACT = "generated_bundle_v2"
CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION = 2
INDEL_GENERATED_BUNDLE_CONTRACT = "generated_bundle_v3"
INDEL_GENERATED_BUNDLE_SCHEMA_VERSION = 3
SUPPORTED_GENERATED_BUNDLE_CONTRACTS = frozenset(
    {
        "generated_bundle_v1",
        CURRENT_GENERATED_BUNDLE_CONTRACT,
        INDEL_GENERATED_BUNDLE_CONTRACT,
    }
)
FIXED_WIDTH_BINNING_ALGORITHM = "fixed_width_decimal_v1"
CANONICAL_CONTIG_CONVENTION = "canonical_without_chr_prefix"
VALIDATION_ARTIFACT_FILENAMES = MappingProxyType(
    {"report": "validation.json", "marker": "validation.ok"}
)
BUNDLE_COMPONENT_FILENAMES: Mapping[str, str] = MappingProxyType(
    {
        "lookup_het": "lookup_het.tsv.gz",
        "lookup_hom": "lookup_hom.tsv.gz",
        "lookup_hem_chrX": "lookup_hem_chrX.tsv.gz",
        "eligible_bed": "eligible_testing_units.bed.gz",
        "calibration_regions": "calibration_regions.bed.gz",
        "excluded_alleles": "excluded_alleles.tsv.gz",
        "candidate_audit": "candidate_audit.tsv.gz",
        "population_af_zero_imputations": "population_af_zero_imputations.tsv.gz",
        "excluded_units": "excluded_units.tsv.gz",
        "metadata": "metadata.json",
        "validation": VALIDATION_ARTIFACT_FILENAMES["report"],
        "marker": VALIDATION_ARTIFACT_FILENAMES["marker"],
    }
)
CORE_COMPONENTS = frozenset(
    {
        "lookup_het",
        "lookup_hom",
        "eligible_testing_units",
        "calibration_regions",
        "excluded_alleles",
        "excluded_units",
    }
)
HEM_COMPONENT = "lookup_hem_chrX"
GENERATED_ROLES = ("het", "hom", "hem_chrX")


def model_roles(hemizygote: bool) -> tuple[str, ...]:
    """Return enabled lookup roles in the published contract order."""
    return GENERATED_ROLES if hemizygote else GENERATED_ROLES[:2]


class BundleContractError(ValueError):
    """Raised when a bundle cannot be read or upgraded under this release."""


def fixed_width_bin_count(score_max: float, bin_width: float) -> int:
    """Return the bin count for the ``fixed_width_decimal_v1`` grid.

    The contract accepts ordinary floating-point metadata, but treats a ratio
    within 64 ulps of an integer boundary as that integer.  This prevents
    binary representations such as ``0.3 / 0.1`` from spuriously allocating
    an extra bin, and is deliberately identical to R's
    ``fixed_width_bin_count`` implementation.
    """
    try:
        maximum = float(score_max)
        width = float(bin_width)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "score_max and bin_width must be positive finite numbers"
        ) from exc
    if (
        not math.isfinite(maximum)
        or maximum <= 0
        or not math.isfinite(width)
        or width <= 0
    ):
        raise ValueError("score_max and bin_width must be positive finite numbers")
    scaled = maximum / width
    if not math.isfinite(scaled):
        raise ValueError("score_max / bin_width must be finite")
    nearest = round(scaled)
    tolerance = 64 * sys.float_info.epsilon * max(1.0, abs(scaled))
    if abs(scaled - nearest) <= tolerance:
        scaled = float(nearest)
    return math.ceil(scaled)


def is_supported_contract(contract: str) -> bool:
    """Return whether a certificate contract has an explicit reader policy."""
    return str(contract) in SUPPORTED_GENERATED_BUNDLE_CONTRACTS


def require_supported_contract(contract: str) -> str:
    """Return a supported contract or raise a migration-focused error."""
    contract = str(contract)
    if not is_supported_contract(contract):
        supported = ", ".join(sorted(SUPPORTED_GENERATED_BUNDLE_CONTRACTS))
        raise BundleContractError(
            f"Unsupported generated bundle contract {contract!r}; supported contracts: {supported}"
        )
    return contract
