"""Gate a workflow branch on a JSON validation report."""

import json
import sys
from pathlib import Path


def main() -> None:
    report_path = Path(snakemake.input.report)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read validation report {report_path}: {exc}") from exc
    if not isinstance(report, dict) or not isinstance(report.get("errors"), list):
        raise RuntimeError(f"Validation report has no errors list: {report_path}")

    errors = report["errors"]
    if errors:
        if snakemake.params.include_errors:
            raise RuntimeError(
                f"{snakemake.params.failure_message}:\n" + "\n".join(errors)
            )
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise RuntimeError(f"{snakemake.params.failure_message}; see {report_path}")

    Path(snakemake.output.ok).write_text(
        snakemake.params.success_marker,
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
