"""Schema validation shared by Snakemake and direct workflow-plan callers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import jsonschema
import yaml


class ConfigurationValidationError(ValueError):
    """Raised when a merged configuration fails its declarative contract."""


SCHEMA_DIR = Path(__file__).resolve().parents[2] / "schemas"
SCHEMA_PATH = SCHEMA_DIR / "config.schema.yaml"


def _load_schema(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def validate_configuration(config: Mapping[str, Any]) -> None:
    """Validate one merged configuration using only repository-local schemas."""
    schema_paths = {
        path.name: path
        for path in (
            SCHEMA_PATH,
            SCHEMA_DIR / "run.schema.yaml",
            SCHEMA_DIR / "science.schema.yaml",
        )
    }
    schemas = {name: _load_schema(path) for name, path in schema_paths.items()}
    schema = schemas[SCHEMA_PATH.name]
    resolver = jsonschema.RefResolver(
        base_uri=SCHEMA_PATH.resolve().as_uri(),
        referrer=schema,
        store={path.resolve().as_uri(): schemas[name] for name, path in schema_paths.items()},
    )
    validator = jsonschema.Draft202012Validator(schema, resolver=resolver)
    def leaf_errors(error: jsonschema.ValidationError) -> list[jsonschema.ValidationError]:
        if not error.context:
            return [error]
        return [leaf for child in error.context for leaf in leaf_errors(child)]

    errors = [
        leaf
        for error in validator.iter_errors(dict(config))
        for leaf in leaf_errors(error)
    ]
    if not errors:
        return

    error = sorted(
        errors,
        key=lambda candidate: (
            candidate.validator != "type",
            -len(candidate.absolute_path),
            list(candidate.absolute_path),
        ),
    )[0]
    path = ".".join(str(part) for part in error.absolute_path) or "<root>"
    if error.validator == "type" and error.validator_value == "boolean":
        raise ConfigurationValidationError(f"{path} must be a boolean")
    raise ConfigurationValidationError(
        f"Configuration schema validation failed at {path}: {error.message}"
    )
