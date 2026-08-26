#!/usr/bin/env python3
"""Merge independent validation fragments into the public validation report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--fragment", action="append", required=True)
    result.add_argument("--output", required=True)
    return result


def append_unique(target: list[str], values: list[str]) -> None:
    seen = set(target)
    for value in values:
        if value not in seen:
            target.append(value)
            seen.add(value)


def main() -> int:
    args = parser().parse_args()
    errors: list[str] = []
    warnings: list[str] = []
    summary: dict = {}
    sections: list[str] = []

    for fragment_text in args.fragment:
        fragment_path = Path(fragment_text)
        fragment = json.loads(fragment_path.read_text(encoding="utf-8"))
        append_unique(errors, fragment.get("errors", []))
        append_unique(warnings, fragment.get("warnings", []))
        summary.update(fragment.get("summary", {}))
        sections.append(fragment_path.stem)

    report = {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": warnings,
        "summary": summary,
        "validation_sections": sections,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"Merged {len(sections)} validation sections into {output}: "
        f"{len(errors)} error(s), {len(warnings)} warning(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
