"""Model and lookup-bundle validation for the apply barrier."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping

from bundle_contract import is_supported_contract
from bundle_certificate import BundleCertificateError, verify_bundle_certificate
from validation import (
    validate_lookup_files,
    validate_lookup_unit_overlap,
)
from validation_reporting import build_validation_report


VALID_UNITS = frozenset({"gene", "cadd_region", "coding_cadd_region"})


def certificate_components_with_exclusions(
    components: Mapping[str, str | Path], excluded_alleles: Iterable[str | Path]
) -> dict[str, str]:
    """Name every declared exclusion input so certificate verification sees it.

    Bundle metadata currently certifies one ``excluded_alleles`` component.
    Additional paths are deliberately given unique component names: the
    certificate verifier will then reject them unless a future contract
    explicitly certifies each one, rather than silently trusting them.
    """
    resolved = {name: str(path) for name, path in components.items()}
    for index, path in enumerate(excluded_alleles):
        name = "excluded_alleles" if index == 0 else f"excluded_alleles_{index:03d}"
        resolved[name] = str(path)
    return resolved


def validate_model_inputs(
    *,
    score_max: float,
    assembly: str,
    unit: str,
    hemizygote_model: bool,
    include_indels: bool,
    testing_unit_bed: str | Path,
    lookup_tables: Mapping[str, str],
    lookup_contract: str,
    lookup_bin_width: float | None,
    bundle_certificate: Mapping[str, str | Path] | None,
    bundle_components: Mapping[str, str | Path],
) -> dict[str, object]:
    """Validate model settings and either an external or certified lookup source."""
    errors: list[str] = []
    warnings: list[str] = []
    if unit not in VALID_UNITS:
        errors.append(
            "unit must be one of "
            + ", ".join(sorted(VALID_UNITS))
            + f"; received {unit!r}"
        )
    if not assembly.strip():
        errors.append("assembly cannot be empty")
    if score_max <= 0:
        errors.append(f"score_max must be positive; received {score_max!r}")
    settings = {
        "hemizygote": hemizygote_model,
        "include_indels": include_indels,
    }
    checked = validate_lookup_files(lookup_tables, settings, errors)
    overlap = validate_lookup_unit_overlap(
        testing_unit_bed,
        lookup_tables,
        settings,
        errors,
    )
    if is_supported_contract(lookup_contract):
        try:
            if bundle_certificate is None:
                raise BundleCertificateError(
                    "Generated lookup bundle is missing its validation certificate"
                )
            metadata = verify_bundle_certificate(
                contract=lookup_contract,
                metadata_path=bundle_certificate["metadata"],
                validation_path=bundle_certificate["validation"],
                marker_path=bundle_certificate["marker"],
                components=bundle_components,
            )
            roles = metadata["model_assumptions"]["roles"]
            variant_classes = metadata["model_assumptions"].get("variant_classes", [])
            bundle_has_indels = "indel" in variant_classes
            if bool(settings["include_indels"]) != bundle_has_indels:
                errors.append(
                    "Generated lookup variant classes do not match model.include_indels"
                )
            if bool(settings["hemizygote"]) != ("hem_chrX" in roles):
                errors.append(
                    "Generated lookup HEM role does not match model.hemizygote"
                )
            if metadata.get("assembly") != assembly:
                errors.append(
                    "Generated lookup assembly differs from reference.assembly"
                )
            if metadata.get("testing_unit_type") != unit:
                errors.append(
                    "Generated lookup testing-unit type differs from model.unit"
                )
            binning = metadata["binning"]
            if not math.isclose(
                float(binning["score_max"]),
                score_max,
                rel_tol=0,
                abs_tol=1e-12,
            ):
                errors.append("Generated lookup score_max differs from model.score_max")
            if lookup_bin_width is None:
                errors.append(
                    "Generated lookup bundle has no resolved lookup bin_width"
                )
            elif not math.isclose(
                float(binning["bin_width"]),
                float(lookup_bin_width),
                rel_tol=0,
                abs_tol=1e-15,
            ):
                errors.append(
                    "Generated lookup bin_width differs from the resolved lookup source"
                )
        except (BundleCertificateError, KeyError, TypeError, ValueError) as exc:
            errors.append(f"Could not validate generated lookup bundle contract: {exc}")
    elif lookup_contract != "external_legacy_v0":
        errors.append(f"Unsupported lookup contract: {lookup_contract}")
    return build_validation_report(
        errors=errors,
        warnings=warnings,
        summary={
            "unit": unit,
            "assembly": assembly,
            "hemizygote_model": settings["hemizygote"],
            "include_indels": settings["include_indels"],
            "lookup_contract": lookup_contract,
            "lookup_score_max": score_max,
            "lookup_bin_width": lookup_bin_width,
            "lookup_files_checked": checked,
            "lookup_testing_unit_overlap": overlap,
        },
        sections=("model",),
    )
