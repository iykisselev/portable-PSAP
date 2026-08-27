"""Write annotation completeness QC using injected Snakemake inputs and outputs."""

import json
from pathlib import Path

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from annotation_qc import summarize_annotations


def main() -> None:
    """Write annotation-completeness QC for the configured cohort annotations.

    Reads the flattened annotation table and indel setting from Snakemake, then
    writes a deterministic JSON report consumed by the annotation validation
    gate. Input-format errors propagate from `summarize_annotations`.
    """
    report = summarize_annotations(
        snakemake.input.annotations,
        include_indels=bool(snakemake.params.include_indels),
    )
    Path(snakemake.output.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
