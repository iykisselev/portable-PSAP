"""Summarize annotation completeness for the apply workflow."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import AnnotationInputs, ReportOutput, AnnotationParameters
from _library import bind


from annotation_qc import summarize_annotations
from validation_reporting import build_validation_report, write_validation_report


def main() -> None:
    """Write annotation-completeness QC for the configured cohort annotations.

    Reads the flattened annotation table and indel setting from Snakemake, then
    writes a deterministic JSON report consumed by the annotation validation
    gate. Input-format errors propagate from `summarize_annotations`.
    """
    inputs = bind(AnnotationInputs, snakemake.input)
    outputs = bind(ReportOutput, snakemake.output)
    params = bind(AnnotationParameters, snakemake.params)
    report = summarize_annotations(
        inputs.annotations,
        include_indels=bool(params.include_indels),
    )
    write_validation_report(
        build_validation_report(
            errors=report["errors"],
            warnings=report["warnings"],
            summary=report["summary"],
            sections=("annotation",),
        ),
        outputs.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
