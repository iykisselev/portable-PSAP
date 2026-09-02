"""Write annotation completeness QC using injected Snakemake inputs and outputs."""

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from annotation_qc import summarize_annotations
from validation_reporting import build_validation_report, write_validation_report


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
    write_validation_report(
        build_validation_report(
            errors=report["errors"],
            warnings=report["warnings"],
            summary=report["summary"],
            sections=("annotation",),
        ),
        snakemake.output.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
