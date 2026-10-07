"""Rendering of bounded statistics into public lookup bucket rows."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)
from generation_manifests import bucket_from_manifest
from generation_probability import distribution, merge_statistics
from io_contracts import write_tsv_atomic


def reduce_partition_statistics(
    inputs: LookupReductionInputs,
    outputs: LookupReductionOutputs,
    params: LookupReductionParameters,
    bucket_id: str,
) -> None:
    """Reduce one unit bucket without materializing genome-wide event tables."""
    units = list(bucket_from_manifest(inputs.buckets, bucket_id).unit_ids)
    grouped = merge_statistics([Path(path) for path in inputs.partition_statistics]) if inputs.partition_statistics else defaultdict(dict)
    lookup_rows = {role: [] for role in params.roles}
    diagnostics = []
    for unit in units:
        for role in params.roles:
            statistics = grouped[(unit, role)]
            aggregate_statistics = {
                bin_number: values[:3] for bin_number, values in statistics.items()
            }
            status, event_mass, event_count, survival = distribution(
                aggregate_statistics, params.bin_count
            )
            snv_event_count = sum(values[3] for values in statistics.values())
            indel_event_count = sum(values[4] for values in statistics.values())
            diagnostics.append(
                (
                    unit,
                    role,
                    event_mass,
                    status,
                    event_count,
                    snv_event_count,
                    indel_event_count,
                )
            )
            if survival is not None:
                lookup_rows[role].append((unit, *survival))
    header = ["testing_unit_id", *[str(value) for value in range(1, params.bin_count + 1)]]
    write_tsv_atomic(outputs.lookup_het, header, lookup_rows["het"])
    write_tsv_atomic(outputs.lookup_hom, header, lookup_rows["hom"])
    if outputs.lookup_hem:
        write_tsv_atomic(outputs.lookup_hem, header, lookup_rows.get("hem_chrX", []))
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
