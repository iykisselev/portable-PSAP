"""The sole public component layout for partitioned generated null bundles."""

from __future__ import annotations

import re
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from model_catalog import enabled_models


PARTITIONED_NULL_BUNDLE_CONTRACT = "partitioned_null_bundle_v1"
PARTITIONED_NULL_BUNDLE_SCHEMA_VERSION = 1
PARTITION_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*\Z")


def require_partition_id(partition_id: str) -> str:
    """Return a safe partition ID for a flat public component filename."""
    if not isinstance(partition_id, str) or not PARTITION_ID_PATTERN.fullmatch(
        partition_id
    ):
        raise ValueError("Partition ID must match ^[A-Za-z0-9][A-Za-z0-9_-]*$")
    return partition_id


def lookup_filename(partition_id: str, model_id: str) -> str:
    require_partition_id(partition_id)
    known_models = {model.identifier for model in enabled_models(True)}
    if model_id not in known_models:
        raise ValueError(f"Unknown public disease model {model_id!r}")
    return f"lookup.{partition_id}.{model_id}.tsv.gz"


def lookup_components(
    partition_ids: tuple[str, ...], hemizygote: bool
) -> Mapping[tuple[str, str], str]:
    """Return every physical partition/model lookup path in config order."""
    if not partition_ids or len(partition_ids) != len(set(partition_ids)):
        raise ValueError("Bundle requires one or more unique partition IDs")
    return MappingProxyType(
        {
            (partition_id, model.identifier): lookup_filename(
                partition_id, model.identifier
            )
            for partition_id in partition_ids
            for model in enabled_models(hemizygote)
        }
    )


def bundle_component_paths(
    directory: str | Path, partition_ids: tuple[str, ...], hemizygote: bool
) -> Mapping[tuple[str, str], str]:
    """Join the public lookup layout to a generated bundle directory."""
    root = Path(directory)
    return MappingProxyType(
        {
            key: str(root / filename)
            for key, filename in lookup_components(partition_ids, hemizygote).items()
        }
    )
