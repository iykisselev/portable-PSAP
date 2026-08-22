"""Merge validation fragments via injected Snakemake I/O."""

import json
from pathlib import Path


def append_unique(target: list[str], values: list[str]) -> None:
    seen = set(target)
    for value in values:
        if value not in seen:
            target.append(value)
            seen.add(value)


def main() -> None:
    errors: list[str] = []
    warnings: list[str] = []
    summary: dict = {}
    sections: list[str] = []
    for fragment_text in snakemake.input:
        fragment_path = Path(fragment_text)
        fragment = json.loads(fragment_path.read_text(encoding="utf-8"))
        append_unique(errors, fragment.get("errors", []))
        append_unique(warnings, fragment.get("warnings", []))
        summary.update(fragment.get("summary", {}))
        sections.append(fragment_path.stem)
    report = {"status": "failed" if errors else "ok", "errors": errors, "warnings": warnings, "summary": summary, "validation_sections": sections}
    Path(snakemake.output.report).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
