"""Create a gate marker only when a validation report has no errors."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import GateInputs, GateOutputs, GateParameters
from _library import bind


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
    inputs = bind(GateInputs, snakemake.input)
    outputs = bind(GateOutputs, snakemake.output)
    params = bind(GateParameters, snakemake.params)
    require_validation_success(
        inputs.report,
        outputs.ok,
        failure_message=params.failure_message,
        include_errors=bool(params.include_errors),
        success_marker=params.success_marker,
    )


if __name__ == "__main__":
    main()
