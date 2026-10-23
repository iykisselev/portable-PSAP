"""Authoritative disease-model mappings and status domains shared with R."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from genomic_contracts import AUTOSOMES, canonical_chromosome


CATALOG_PATH = Path(__file__).resolve().parents[2] / "model_catalog.json"


@dataclass(frozen=True)
class DiseaseModel:
    identifier: str
    role: str
    scope: str
    genotype_class: str
    tag: str
    requires_hemizygote: bool
    compatible_control_models: tuple[str, ...]


@dataclass(frozen=True)
class StatusDefinition:
    code: str
    domains: tuple[str, ...]
    required_payload: Mapping[str, tuple[str, ...]]
    diagnostic_only: bool


@dataclass(frozen=True)
class ModelCatalog:
    roles: tuple[str, ...]
    models: tuple[DiseaseModel, ...]
    statuses: Mapping[str, StatusDefinition]


def _strings(value: object, label: str, *, nonempty: bool = True) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or (nonempty and not value)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"Catalog {label} must be an array of unique nonempty strings")
    return tuple(value)


def read_model_catalog(path: str | Path = CATALOG_PATH) -> ModelCatalog:
    """Read a validated immutable catalog; never supply fallback model mappings."""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ValueError("Unsupported model catalog schema")
    roles = _strings(value.get("roles"), "roles")
    rows = value.get("models")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Catalog models must be a nonempty array")
    models: list[DiseaseModel] = []
    identifiers: set[str] = set()
    tags: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Catalog model entries must be objects")
        for name in ("id", "role", "scope", "genotype_class", "tag"):
            if not isinstance(row.get(name), str) or not row[name]:
                raise ValueError(f"Catalog model requires {name}")
        if (
            row["id"] in identifiers
            or row["tag"] in tags
            or row["role"] not in roles
            or row["scope"] not in {"autosomal", "non_PAR_chrX"}
            or not isinstance(row.get("requires_hemizygote"), bool)
        ):
            raise ValueError(f"Invalid or duplicate catalog model {row['id']!r}")
        identifiers.add(row["id"])
        tags.add(row["tag"])
        models.append(
            DiseaseModel(
                identifier=row["id"],
                role=row["role"],
                scope=row["scope"],
                genotype_class=row["genotype_class"],
                tag=row["tag"],
                requires_hemizygote=row["requires_hemizygote"],
                compatible_control_models=_strings(
                    row.get("compatible_control_models"), "compatible_control_models"
                ),
            )
        )
    if any(set(model.compatible_control_models) - identifiers for model in models):
        raise ValueError("Catalog contains an unknown compatible control model")
    status_rows = value.get("statuses")
    if not isinstance(status_rows, list) or not status_rows:
        raise ValueError("Catalog statuses must be a nonempty array")
    statuses: dict[str, StatusDefinition] = {}
    for row in status_rows:
        if not isinstance(row, dict):
            raise ValueError("Catalog status entries must be objects")
        code = row.get("code")
        if not isinstance(code, str) or not code or code in statuses:
            raise ValueError("Invalid or duplicate catalog status")
        domains = _strings(row.get("domains"), "status domains")
        if set(domains) - {"scope", "role", "partition", "aggregation", "triage"}:
            raise ValueError(f"Unknown catalog status domain for {code}")
        required = row.get("required_payload")
        if not isinstance(required, dict) or set(required) != set(domains):
            raise ValueError(
                f"Catalog status {code} must declare payloads for every domain"
            )
        payloads = {
            domain: _strings(fields, "required_payload", nonempty=False)
            for domain, fields in required.items()
        }
        if any(set(fields) - {"role", "failed_roles"} for fields in payloads.values()):
            raise ValueError(f"Unknown catalog payload field for {code}")
        if not isinstance(row.get("diagnostic_only"), bool):
            raise ValueError(f"Catalog status {code} requires diagnostic_only boolean")
        statuses[code] = StatusDefinition(
            code, domains, MappingProxyType(payloads), row["diagnostic_only"]
        )
    return ModelCatalog(roles, tuple(models), MappingProxyType(statuses))


CATALOG = read_model_catalog()
GENERATED_ROLES = CATALOG.roles


def enabled_models(hemizygote: bool) -> tuple[DiseaseModel, ...]:
    if not isinstance(hemizygote, bool):
        raise ValueError("hemizygote must be a boolean")
    return tuple(m for m in CATALOG.models if hemizygote or not m.requires_hemizygote)


def model_roles(hemizygote: bool) -> tuple[str, ...]:
    """Return roles in catalog order, selecting them from the enabled models."""
    used = {model.role for model in enabled_models(hemizygote)}
    return tuple(role for role in CATALOG.roles if role in used)


def models_for_chromosome(
    chromosome: str, hemizygote: bool
) -> tuple[DiseaseModel, ...]:
    chromosome = canonical_chromosome(chromosome)
    scope = (
        "autosomal"
        if chromosome in AUTOSOMES
        else "non_PAR_chrX"
        if chromosome == "X"
        else None
    )
    return tuple(model for model in enabled_models(hemizygote) if model.scope == scope)


def validate_status(
    code: str,
    domain: str,
    *,
    role: str | None = None,
    failed_roles: tuple[str, ...] = (),
    production: bool = True,
) -> str:
    """Validate a status in its domain, including structured role diagnostics."""
    definition = CATALOG.statuses.get(code)
    if definition is None or domain not in definition.domains:
        raise ValueError(f"Unsupported {domain} status {code!r}")
    if production and definition.diagnostic_only:
        raise ValueError(f"Status {code!r} is diagnostic only")
    if role is not None and role not in CATALOG.roles:
        raise ValueError(f"Unknown diagnostic role {role!r}")
    if (
        not isinstance(failed_roles, (tuple, list))
        or any(item not in CATALOG.roles for item in failed_roles)
        or len(failed_roles) != len(set(failed_roles))
        or tuple(failed_roles) != tuple(r for r in CATALOG.roles if r in failed_roles)
    ):
        raise ValueError("failed_roles must contain unique roles in catalog order")
    required = definition.required_payload[domain]
    if "role" in required and role is None:
        raise ValueError(f"Status {code!r} requires role")
    if "failed_roles" in required and not failed_roles:
        raise ValueError(f"Status {code!r} requires nonempty failed_roles")
    if failed_roles and "failed_roles" not in required:
        raise ValueError(f"Status {code!r} does not allow failed_roles in {domain}")
    if role is not None and "role" not in required:
        raise ValueError(f"Status {code!r} does not allow role in {domain}")
    return code
