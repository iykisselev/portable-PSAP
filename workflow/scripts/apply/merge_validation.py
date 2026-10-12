"""Merge validation fragments into the single report used by apply rules."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import MergeInputs, ReportOutput
from _library import bind


from validation_reporting import merge_validation_reports, write_validation_report


def main() -> None:
    """Merge validation fragments into the workflow's single gate report.

    Reads JSON reports listed by Snakemake and writes one report containing
    ordered, de-duplicated errors and warnings. Summary keys follow input order,
    so a later fragment intentionally replaces an earlier value with the same
    key.

    Raises:
        OSError: If a fragment or the destination cannot be read or written.
        json.JSONDecodeError: If a validation fragment is not valid JSON.
    """
    inputs = bind(MergeInputs, snakemake.input)
    outputs = bind(ReportOutput, snakemake.output)
    write_validation_report(
        merge_validation_reports((inputs.cohort, inputs.resources, inputs.model)),
        outputs.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
