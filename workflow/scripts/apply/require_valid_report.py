"""Create a gate marker only when a validation report has no errors."""

from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from validation_reporting import require_validation_success


def main() -> None:
    """Create a gate marker only when a validation report has no errors.

    Reads the JSON report supplied by Snakemake. On success it writes the
    configured marker; on failure it emits errors to stderr or includes them in
    the raised message according to the rule configuration.

    Raises:
        RuntimeError: If the report cannot be parsed, lacks an errors list, or
            records validation failures.
    """
    require_validation_success(
        snakemake.input.report,
        snakemake.output.ok,
        failure_message=snakemake.params.failure_message,
        include_errors=bool(snakemake.params.include_errors),
        success_marker=snakemake.params.success_marker,
    )


if __name__ == "__main__":
    main()
