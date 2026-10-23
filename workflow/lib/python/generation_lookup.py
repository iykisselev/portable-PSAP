"""Rendering of bounded statistics into public lookup bucket rows."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)
from generation_manifests import bucket_from_manifest
from generation_probability import aggregate_statistics_vectorized
from io_contracts import write_tsv_atomic
from partitioned_lookup import partition_model_lookup_rows


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
    header = [
        "testing_unit_id",
        *[str(value) for value in range(1, params.bin_count + 1)],
    ]
    if len(outputs.partition_model_lookups) != len(params.partition_lookup_keys):
        raise ValueError("Partition/model bucket lookup outputs do not match keys")
    if inputs.partition_statuses is None:
        raise ValueError("Partition/model lookup reduction requires partition statuses")
    pairs = (
        pl.scan_parquet(inputs.partition_statuses)
        .filter(
            (pl.col("status") == "eligible") & pl.col("testing_unit_id").is_in(units)
        )
        .select("testing_unit_id", "partition_id", "CHROM")
        .collect(engine="streaming")
        .to_dicts()
    )
    model_rows, model_diagnostics = partition_model_lookup_rows(
        statistics,
        pairs,
        hemizygote="hem_chrX" in params.roles,
        bin_count=params.bin_count,
    )
    for path, key in zip(
        outputs.partition_model_lookups, params.partition_lookup_keys, strict=True
    ):
        partition_id, model_id = key.split("\t", 1)
        write_tsv_atomic(path, header, model_rows[(partition_id, model_id)])
    write_tsv_atomic(
        outputs.partition_model_diagnostics,
        [
            "partition_id",
            "model",
            "testing_unit_id",
            "role",
            "event_mass",
            "status",
            "event_count",
            "snv_event_count",
            "indel_event_count",
        ],
        (
            [
                row[column]
                for column in (
                    "partition_id",
                    "model",
                    "testing_unit_id",
                    "role",
                    "event_mass",
                    "status",
                    "event_count",
                    "snv_event_count",
                    "indel_event_count",
                )
            ]
            for row in model_diagnostics
        ),
    )
