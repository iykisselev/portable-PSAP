"""Write annotation completeness QC using injected Snakemake inputs and outputs."""

import json
from pathlib import Path

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from annotation_qc import summarize_annotations


def main() -> None:
    report = summarize_annotations(
        snakemake.input.annotations,
        include_indels=bool(snakemake.params.include_indels),
    )
    Path(snakemake.output.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
