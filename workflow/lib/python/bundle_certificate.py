"""Shared verification of generated null-distribution bundle certificates."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

try:
    from bundle_contract import require_supported_contract
    from io_contracts import sha256_file
except ModuleNotFoundError:
    from .bundle_contract import require_supported_contract
    from .io_contracts import sha256_file


class BundleCertificateError(ValueError):
    """Raised when a bundle's published certificate cannot be trusted."""


def validate_bundle_component(
    name: str,
    path: str | Path,
    metadata: Mapping[str, object],
) -> None:
    """Check the filename, byte count, and digest declared for one component."""
    try:
        declared = metadata["components"][name]
        component = Path(path)
        if component.name != declared["path"]:
            raise BundleCertificateError(
                f"Generated bundle component {name} has the wrong filename"
            )
        if component.stat().st_size != declared["bytes"]:
            raise BundleCertificateError(
                f"Generated bundle component {name} has the wrong size"
            )
        if sha256_file(component) != declared["sha256"]:
            raise BundleCertificateError(
                f"Generated bundle component {name} checksum does not match metadata"
            )
    except (KeyError, OSError, TypeError) as exc:
        raise BundleCertificateError(
            f"Could not validate generated bundle component {name}: {exc}"
        ) from exc


def verify_bundle_certificate(
    *,
    contract: str,
    metadata_path: str | Path,
    validation_path: str | Path,
    marker_path: str | Path,
    components: Mapping[str, str | Path],
) -> dict[str, object]:
    """Verify one producer-issued certificate and each bound public component."""
    require_supported_contract(contract)
    try:
        metadata_file = Path(metadata_path)
        validation_file = Path(validation_path)
        marker_file = Path(marker_path)
        metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
        validation = json.loads(validation_file.read_text(encoding="utf-8"))
        metadata_digest = sha256_file(metadata_file)
        validation_digest = sha256_file(validation_file)
        marker_fields = marker_file.read_text(encoding="utf-8").strip().split("\t")
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise BundleCertificateError(
            f"Could not read generated lookup bundle certificate: {exc}"
        ) from exc

    if metadata.get("contract") != contract:
        raise BundleCertificateError("Generated lookup metadata has the wrong contract")
    if validation.get("status") != "ok" or validation.get("metadata_sha256") != metadata_digest:
        raise BundleCertificateError(
            "Generated lookup validation does not certify metadata.json"
        )
    if marker_fields != [contract, metadata_digest, validation_digest]:
        raise BundleCertificateError(
            "Generated lookup marker does not certify metadata.json and validation.json"
        )
    for name, path in components.items():
        validate_bundle_component(name, path, metadata)
    return metadata
