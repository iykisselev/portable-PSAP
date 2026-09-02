"""Merge schema-validated fragments into the one apply validation report."""

from _library import add_library_path

add_library_path(snakemake.scriptdir)

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
    write_validation_report(
        merge_validation_reports(snakemake.input),
        snakemake.output.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
