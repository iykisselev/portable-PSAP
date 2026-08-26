#!/usr/bin/env python3
"""Fail the workflow when the persisted input-validation report contains errors."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Require a successful PSAP input-validation report."
    )
    parser.add_argument("--report", required=True)
    parser.add_argument("--output", required=True)
    return parser


def load_report(path: str | Path) -> dict:
    report_path = Path(path)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read validation report {report_path}: {exc}") from exc

    if not isinstance(report, dict):
        raise ValueError(f"Validation report must contain a JSON object: {report_path}")
    if not isinstance(report.get("errors"), list):
        raise ValueError(f"Validation report has no errors list: {report_path}")

    return report


def main() -> int:
    args = build_parser().parse_args()

    try:
        report = load_report(args.report)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if report["errors"]:
        print(
            f"Input validation failed with {len(report['errors'])} error(s).",
            file=sys.stderr,
        )
        for error in report["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"Full validation report: {args.report}", file=sys.stderr)
        return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("validated\n", encoding="utf-8")
    print(f"Input validation gate passed: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
