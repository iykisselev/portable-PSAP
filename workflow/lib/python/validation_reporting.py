"""Schema-stable, atomic reporting primitives for every validation barrier."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Iterable, Mapping

from io_contracts import write_text_atomic


VALIDATION_REPORT_SCHEMA_VERSION = 1


class ValidationReportError(ValueError):
    """Raised when a validation report cannot safely drive a workflow gate."""


def build_validation_report(
    *,
    errors: Iterable[str],
    warnings: Iterable[str],
    summary: Mapping[str, object],
    sections: Iterable[str] = (),
) -> dict[str, object]:
    """Construct one report conforming to the shared validation contract."""
    collected_errors = [str(error) for error in errors]
    return {
        "schema_version": VALIDATION_REPORT_SCHEMA_VERSION,
        "status": "failed" if collected_errors else "ok",
        "errors": collected_errors,
        "warnings": [str(warning) for warning in warnings],
        "summary": dict(summary),
        "validation_sections": [str(section) for section in sections],
    }


def validate_validation_report(report: object) -> dict[str, object]:
    """Validate the small report shape consumed by marker and merge rules."""
    if not isinstance(report, dict):
        raise ValidationReportError("Validation report must be a JSON object")
    required = {
        "schema_version",
        "status",
        "errors",
        "warnings",
        "summary",
        "validation_sections",
    }
    missing = sorted(required - set(report))
    if missing:
        raise ValidationReportError(
            "Validation report is missing required field(s): " + ", ".join(missing)
        )
    if report["schema_version"] != VALIDATION_REPORT_SCHEMA_VERSION:
        raise ValidationReportError(
            "Unsupported validation report schema version: "
            f"{report['schema_version']!r}"
        )
    if report["status"] not in {"ok", "failed"}:
        raise ValidationReportError("Validation report status must be 'ok' or 'failed'")
    for name in ("errors", "warnings", "validation_sections"):
        if not isinstance(report[name], list) or not all(
            isinstance(value, str) for value in report[name]
        ):
            raise ValidationReportError(
                f"Validation report {name} must be a list of strings"
            )
    if not isinstance(report["summary"], dict):
        raise ValidationReportError("Validation report summary must be an object")
    expected_status = "failed" if report["errors"] else "ok"
    if report["status"] != expected_status:
        raise ValidationReportError(
            "Validation report status does not agree with its errors list"
        )
    return report


def write_validation_report(
    report: Mapping[str, object], path: str | Path, *, create_parents: bool = True
) -> Path:
    """Validate and atomically publish a report consumed by a workflow marker."""
    checked = validate_validation_report(dict(report))
    return write_text_atomic(
        path,
        json.dumps(checked, indent=2, sort_keys=True) + "\n",
        create_parents=create_parents,
    )


def read_validation_report(path: str | Path) -> dict[str, object]:
    """Read and validate one report before merging or creating a success marker."""
    source = Path(path)
    try:
        report = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationReportError(
            f"Could not read validation report {source}: {exc}"
        ) from exc
    return validate_validation_report(report)


def merge_validation_reports(paths: Iterable[str | Path]) -> dict[str, object]:
    """Merge independent category reports with ordered, unique diagnostics."""
    errors: list[str] = []
    warnings: list[str] = []
    summary: dict[str, object] = {}
    sections: list[str] = []
    for path in paths:
        report = read_validation_report(path)
        _append_unique(errors, report["errors"])
        _append_unique(warnings, report["warnings"])
        summary.update(report["summary"])
        source_sections = report["validation_sections"] or [Path(path).stem]
        _append_unique(sections, source_sections)
    return build_validation_report(
        errors=errors,
        warnings=warnings,
        summary=summary,
        sections=sections,
    )


def require_validation_success(
    report_path: str | Path,
    marker_path: str | Path,
    *,
    failure_message: str,
    include_errors: bool,
    success_marker: str,
) -> Path:
    """Fail with report diagnostics or atomically write a named success marker."""
    report = read_validation_report(report_path)
    errors = report["errors"]
    if errors:
        if include_errors:
            raise RuntimeError(f"{failure_message}:\n" + "\n".join(errors))
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise RuntimeError(f"{failure_message}; see {report_path}")
    return write_text_atomic(marker_path, success_marker, create_parents=True)


def _append_unique(target: list[str], values: Iterable[object]) -> None:
    seen = set(target)
    for value in values:
        message = str(value)
        if message not in seen:
            target.append(message)
            seen.add(message)
