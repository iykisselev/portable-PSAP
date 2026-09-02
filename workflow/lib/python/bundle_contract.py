"""Versioned generated-bundle names, roles, and upgrade policy.

Only :data:`CURRENT_GENERATED_BUNDLE_CONTRACT` may be emitted by this version
of the workflow.  Older contracts are accepted by consumers through the
explicit compatibility registry; a future writer must add an upgrader rather
than silently reinterpret an existing certificate.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Callable, Mapping


CURRENT_GENERATED_BUNDLE_CONTRACT = "generated_bundle_v2"
CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION = 2
SUPPORTED_GENERATED_BUNDLE_CONTRACTS = frozenset(
    {"generated_bundle_v1", CURRENT_GENERATED_BUNDLE_CONTRACT}
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
        "lookup_hem": "lookup_hem_chrX.tsv.gz",
        "eligible_bed": "eligible_testing_units.bed.gz",
        "calibration_regions": "calibration_regions.bed.gz",
        "excluded_alleles": "excluded_alleles.tsv.gz",
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


class BundleContractError(ValueError):
    """Raised when a bundle cannot be read or upgraded under this release."""


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


def upgrade_metadata(metadata: Mapping[str, object], target: str = CURRENT_GENERATED_BUNDLE_CONTRACT) -> dict[str, object]:
    """Upgrade metadata only through a declared conversion path.

    No v1-to-v2 metadata-only conversion exists: v2 adds certified X-policy
    semantics and potentially an HEM component, so users must regenerate a v1
    bundle.  Keeping that decision executable prevents accidental silent
    contract drift when a later version adds an upgrader.
    """
    source = require_supported_contract(str(metadata.get("contract", "")))
    target = str(target)
    if source == target:
        return dict(metadata)
    raise BundleContractError(
        f"No safe metadata-only upgrade path exists from {source} to {target}; regenerate the bundle"
    )
