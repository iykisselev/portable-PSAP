"""Provenance hashing for raw generation inputs."""

from __future__ import annotations

from collections.abc import Mapping

from generation_bindings import GenerationResourceIdentityOutputs, PathLike
from generation_manifests import PARTITION_MANIFEST_SCHEMA_VERSION
from io_contracts import resource_inventory, write_json_atomic


def hash_generation_resources(
    resources: Mapping[str, PathLike],
    outputs: GenerationResourceIdentityOutputs,
) -> None:
    """Publish source paths and SHA-256 digests independently of the checkpoint."""
    write_json_atomic(
        outputs.resource_identities,
        {
            "schema_version": PARTITION_MANIFEST_SCHEMA_VERSION,
            "resources": resource_inventory(resources),
        },
    )
