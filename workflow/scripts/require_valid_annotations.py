"""Fail the workflow when annotation completeness QC reports an error."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    report_path = Path(snakemake.input.report)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    errors = report.get("errors", [])
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise RuntimeError(f"Annotation validation failed; see {report_path}")
    Path(snakemake.output.ok).write_text("annotation validation passed\n", encoding="utf-8")


if __name__ == "__main__":
    main()
