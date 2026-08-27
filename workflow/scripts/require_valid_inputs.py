"""Gate execution on an injected input-validation report."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    report_path = Path(snakemake.input.report)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read validation report {report_path}: {exc}") from exc
    if not isinstance(report, dict) or not isinstance(report.get("errors"), list):
        raise RuntimeError(f"Validation report has no errors list: {report_path}")
    if report["errors"]:
        raise RuntimeError("Input validation failed:\n" + "\n".join(report["errors"]))
    Path(snakemake.output.ok).write_text("validated\n", encoding="utf-8")


if __name__ == "__main__":
    main()
