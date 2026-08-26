#!/usr/bin/env python3
"""Gate scoring on a successful annotation-completeness report."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    report_path = Path(args.report)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    errors = report.get("errors", [])
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            f"Annotation validation failed; see {report_path}",
            file=sys.stderr,
        )
        return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("annotation validation passed\n", encoding="utf-8")
    print(f"Annotation validation passed; wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
