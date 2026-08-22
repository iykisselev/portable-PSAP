"""Write annotation completeness QC using injected Snakemake inputs and outputs."""

import importlib.util
import json
from pathlib import Path


def _implementation():
    path = Path(__file__).parent / "apply" / "summarize_annotations.py"
    spec = importlib.util.spec_from_file_location("annotation_summary", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def main() -> None:
    implementation = _implementation()
    report = implementation.summarize_annotations(
        snakemake.input.annotations,
        include_indels=bool(snakemake.params.include_indels),
    )
    Path(snakemake.output.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
