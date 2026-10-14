"""Rendering of bounded statistics into public lookup bucket rows."""

from __future__ import annotations

from pathlib import Path

from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)
from generation_manifests import bucket_from_manifest
from generation_probability import (
    aggregate_statistics_vectorized,
    vectorized_lookup_rows,
)
from io_contracts import write_tsv_atomic


def reduce_window_statistics(
    inputs: LookupReductionInputs,
    outputs: LookupReductionOutputs,
    params: LookupReductionParameters,
    bucket_id: str,
) -> None:
    """Reduce one bounded unit bucket through Polars and dense NumPy arrays."""
    units = list(bucket_from_manifest(inputs.buckets, bucket_id).unit_ids)
    statistics = aggregate_statistics_vectorized(
        [Path(path) for path in inputs.window_statistics]
    )
    lookup_rows, diagnostics = vectorized_lookup_rows(
        statistics, units, params.roles, params.bin_count
    )
    header = [
        "testing_unit_id",
        *[str(value) for value in range(1, params.bin_count + 1)],
    ]
    write_tsv_atomic(outputs.lookup_het, header, lookup_rows["het"])
    write_tsv_atomic(outputs.lookup_hom, header, lookup_rows["hom"])
    if outputs.lookup_hem_chrX:
        write_tsv_atomic(
            outputs.lookup_hem_chrX, header, lookup_rows.get("hem_chrX", [])
        )
    write_tsv_atomic(
        outputs.diagnostics,
        [
            "testing_unit_id",
            "role",
            "event_mass",
            "status",
            "event_count",
            "snv_event_count",
            "indel_event_count",
        ],
        diagnostics,
    )
