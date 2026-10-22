"""Shared constants and validation for the partitioned-null bundle contract."""

from __future__ import annotations

import math
import sys
from types import MappingProxyType
from typing import Mapping

from model_catalog import model_roles  # noqa: F401


CURRENT_GENERATED_BUNDLE_CONTRACT = "partitioned_null_bundle_v1"
CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION = 1
INDEL_GENERATED_BUNDLE_CONTRACT = CURRENT_GENERATED_BUNDLE_CONTRACT
INDEL_GENERATED_BUNDLE_SCHEMA_VERSION = CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION
SUPPORTED_GENERATED_BUNDLE_CONTRACTS = frozenset({CURRENT_GENERATED_BUNDLE_CONTRACT})
FIXED_WIDTH_BINNING_ALGORITHM = "fixed_width_decimal_v1"
CANONICAL_CONTIG_CONVENTION = "canonical_without_chr_prefix"
VALIDATION_ARTIFACT_FILENAMES = MappingProxyType(
    {"report": "validation.json", "marker": "validation.ok"}
)
BUNDLE_ARTIFACT_FILENAMES: Mapping[str, str] = MappingProxyType(
    {
        "metadata": "metadata.json",
        "validation": VALIDATION_ARTIFACT_FILENAMES["report"],
        "marker": VALIDATION_ARTIFACT_FILENAMES["marker"],
    }
)
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
