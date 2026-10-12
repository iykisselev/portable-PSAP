"""Validate resource contracts through injected Snakemake inputs and parameters."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import (
    ResourceValidationInputs,
    ReportOutput,
    ResourceValidationParameters,
)
from _library import bind, optional_path


from validation import IndexedResource
from validation_reporting import write_validation_report
from validation_resources import validate_resource_inputs


def indexed_resource(namespace, fields) -> IndexedResource | None:
    """Resolve an identity's paired files from the actual job inputs."""
    path = optional_path(getattr(namespace, fields["data"]))
    if path is None:
        return None
    return IndexedResource(path, optional_path(getattr(namespace, fields["index"])))


def main() -> None:
    """Validate reference-resource formats and chromosome compatibility.

    Converts Snakemake bindings into the shared validation contract, checks
    required indexes and annotation resources, and writes the resulting JSON
    report for the workflow validation gate.
    """
    inputs = bind(ResourceValidationInputs, snakemake.input)
    outputs = bind(ReportOutput, snakemake.output)
    params = bind(ResourceValidationParameters, snakemake.params)
    resources = {
        name: indexed_resource(snakemake.input, fields)
        for name, fields in params.resource_fields.items()
    }

    def by_role(role):
        return [
            resources[name]
            for name, fields in params.resource_fields.items()
            if fields["role"] == role and resources[name] is not None
        ]

    write_validation_report(
        validate_resource_inputs(
            vcf=inputs.vcf,
            calibration_regions=resources.get("calibration_regions"),
            exclusions=by_role("variants_exclude"),
            testing_units=resources["testing_units"],
            par_regions=by_role("par_regions"),
            cadd_snv=resources["cadd_snv"],
            cadd_indel=resources.get("cadd_indel"),
            unit=params.unit,
            include_indels=params.include_indels,
        ),
        outputs.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
