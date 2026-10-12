"""Small shared utilities for thin Snakemake script adapters."""

from collections.abc import Mapping
from dataclasses import MISSING, fields
from pathlib import Path
import sys
import types
from typing import Any, Union, get_args, get_origin, get_type_hints


class BindingError(TypeError):
    """Raised when a Snakemake namespace cannot satisfy a typed binding."""


def workflow_root(scriptdir: str | Path) -> Path:
    """Find the workflow root from either a flat or nested script directory."""
    start = Path(scriptdir).resolve()
    for candidate in (start, *start.parents):
        if (candidate / "lib" / "python").is_dir():
            return candidate
    raise ValueError(
        f"Could not find workflow/lib/python above script directory {start}"
    )


def add_library_path(scriptdir: str | Path) -> None:
    """Prepend the workflow's shared Python library to the import path.

    Args:
        scriptdir: Directory containing the executing Snakemake script.

    Side Effects:
        Inserts ``workflow/lib/python`` at the front of ``sys.path`` when it is
        not already present, so local workflow modules take precedence.
    """
    library_dir = workflow_root(scriptdir) / "lib" / "python"
    library = str(library_dir)
    if library not in sys.path:
        sys.path.insert(0, library)


def bind(cls: type, namespace: object):
    """Construct a dataclass from a typed Snakemake input/output/params namespace.

    Required fields absent from the namespace fail here, at the script boundary,
    rather than later inside the generation library.  Dataclass defaults remain
    authoritative for omitted optional rule bindings.
    """
    hints = get_type_hints(cls)
    values: dict[str, object] = {}
    for field in fields(cls):
        if not hasattr(namespace, field.name):
            if field.default is MISSING and field.default_factory is MISSING:
                raise BindingError(
                    f"{cls.__name__} requires Snakemake field {field.name!r}; "
                    f"{type(namespace).__name__} does not provide it"
                )
            continue
        values[field.name] = _coerce(getattr(namespace, field.name), hints[field.name])
    return cls(**values)


def optional_path(value: object) -> str | None:
    """Normalize a Snakemake optional-path value to a string or ``None``."""
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        return str(value) or None
    try:
        if len(value) == 0:  # type: ignore[arg-type]
            return None
        return str(value[0])  # type: ignore[index]
    except TypeError as error:
        raise BindingError(f"Expected an optional path, got {value!r}") from error


def _coerce(value: object, hint: object) -> object:
    origin = get_origin(hint)
    if origin in (Union, types.UnionType):
        options = get_args(hint)
        if type(None) in options:
            path_options = tuple(
                option for option in options if option is not type(None)
            )
            if _is_path_like(path_options):
                return optional_path(value)
            if value is None:
                return None
            return _coerce(value, path_options[0])
        if _is_path_like(options):
            return _path(value)

    if _is_path_like((hint,)):
        return _path(value)
    if origin is tuple:
        item_hint = get_args(hint)[0] if get_args(hint) else Any
        return tuple(_coerce(item, item_hint) for item in value)  # type: ignore[union-attr]
    if origin in (dict, Mapping):
        key_hint, item_hint = get_args(hint) or (Any, Any)
        return {
            _coerce(key, key_hint): _coerce(item, item_hint)
            for key, item in dict(value).items()  # type: ignore[arg-type]
        }
    if hint is str:
        return str(value)
    if hint is int:
        return int(value)
    if hint is float:
        return float(value)
    if hint is bool:
        return _bool(value)
    return value


def _is_path_like(options: tuple[object, ...]) -> bool:
    return bool(options) and set(options) <= {str, Path}


def _path(value: object) -> str:
    if isinstance(value, (str, Path)):
        return str(value)
    try:
        if len(value) == 1:  # type: ignore[arg-type]
            return str(value[0])  # type: ignore[index]
    except TypeError:
        pass
    raise BindingError(f"Expected one path, got {value!r}")


def _bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no"}:
            return False
        raise BindingError(f"Expected a boolean value, got {value!r}")
    return bool(value)
