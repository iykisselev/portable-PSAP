"""Stream per-lookup endpoint diagnostics for the ACAT policy decision."""

import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from _library import bind
from apply_bindings import (
    EndpointDiagnosticsInputs,
    EndpointDiagnosticsOutputs,
    EndpointDiagnosticsParameters,
)
from io_contracts import GzipTsvAtomicWriter, open_text


FIELDS = (
    "partition_id",
    "model",
    "testing_unit_id",
    "bin_count",
    "exact_zero_bin_count",
    "exact_one_bin_count",
    "near_one_0_999_bin_count",
    "near_one_0_999999_bin_count",
    "near_one_0_999999999_bin_count",
    "exact_zero_bin_fraction",
    "exact_one_bin_fraction",
    "near_one_0_999_fraction",
    "near_one_0_999999_fraction",
    "near_one_0_999999999_fraction",
    "exact_zero_carrier_rate",
    "exact_one_carrier_rate",
    "near_one_0_999_carrier_rate",
    "near_one_0_999999_carrier_rate",
    "near_one_0_999999999_carrier_rate",
)


def _row_metrics(
    values: list[float], thresholds: tuple[float, ...]
) -> tuple[object, ...]:
    if not values:
        raise ValueError("Lookup row has no PSAP bins")
    if any(value < 0 or value > 1 for value in values) or any(
        left < right for left, right in zip(values, values[1:], strict=False)
    ):
        raise ValueError("Lookup row has invalid survival values")
    masses = [
        left - right for left, right in zip(values, [*values[1:], 0.0], strict=True)
    ]
    if any(mass < -1e-12 for mass in masses) or abs(sum(masses) - 1.0) > 1e-9:
        raise ValueError("Lookup survival row has invalid event mass")
    # Small arithmetic negatives are only accepted within the mass tolerance.
    masses = [max(0.0, mass) for mass in masses]
    predicates = [
        lambda value: value == 0.0,
        lambda value: value == 1.0,
        *(
            lambda value, threshold=threshold: value >= threshold
            for threshold in thresholds
        ),
    ]
    counts = [sum(predicate(value) for value in values) for predicate in predicates]
    rates = [
        sum(
            mass for mass, value in zip(masses, values, strict=True) if predicate(value)
        )
        for predicate in predicates
    ]
    n = len(values)
    return (*counts, *(count / n for count in counts), *rates)


def main() -> None:
    inputs = bind(EndpointDiagnosticsInputs, snakemake.input)
    outputs = bind(EndpointDiagnosticsOutputs, snakemake.output)
    params = bind(EndpointDiagnosticsParameters, snakemake.params)
    if len(inputs.lookup_tables) != len(params.lookup_keys):
        raise ValueError("Endpoint diagnostic lookup paths and keys do not match")
    if tuple(params.near_one_thresholds) != (0.999, 0.999999, 0.999999999):
        raise ValueError(
            "Endpoint diagnostics require the declared near-one thresholds"
        )
    with GzipTsvAtomicWriter(outputs.diagnostics, FIELDS) as output:
        for key, path in zip(params.lookup_keys, inputs.lookup_tables, strict=True):
            try:
                partition_id, model = key.split("\t", 1)
            except ValueError as exc:
                raise ValueError("Endpoint diagnostic lookup key is invalid") from exc
            with open_text(path, newline="", require_gzip=True) as handle:
                rows = csv.reader(handle, delimiter="\t")
                header = next(rows, None)
                if not header or header[0] != "testing_unit_id":
                    raise ValueError("Endpoint diagnostic lookup has invalid header")
                for row in rows:
                    if len(row) != len(header) or not row[0]:
                        raise ValueError("Endpoint diagnostic lookup has invalid row")
                    try:
                        values = [float(value) for value in row[1:]]
                    except ValueError as exc:
                        raise ValueError(
                            "Endpoint diagnostic lookup has nonnumeric values"
                        ) from exc
                    output.write_row(
                        (
                            partition_id,
                            model,
                            row[0],
                            len(values),
                            *_row_metrics(values, params.near_one_thresholds),
                        )
                    )


if __name__ == "__main__":
    main()
