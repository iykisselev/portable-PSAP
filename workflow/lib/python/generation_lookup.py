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
from generation_probability import (
    aggregate_statistics_vectorized,
    vectorized_lookup_rows,
)
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
    # The old role-wide files cannot represent the same testing unit in more
    # than one biological partition.  Keep them as empty compatibility
    # artifacts only while the Snakemake graph is migrated; publication uses
    # the partition/model rows below.
    if statistics.get_column("partition_id").n_unique() <= 1:
        lookup_rows, diagnostics = vectorized_lookup_rows(
            statistics.drop("partition_id"), units, params.roles, params.bin_count
        )
    else:
        lookup_rows = {role: [] for role in params.roles}
        diagnostics = []
    write_tsv_atomic(outputs.lookup_het, header, lookup_rows["het"])
    write_tsv_atomic(outputs.lookup_hom, header, lookup_rows["hom"])
    if outputs.lookup_hem_chrX:
        write_tsv_atomic(outputs.lookup_hem_chrX, header, lookup_rows.get("hem_chrX", []))
    write_tsv_atomic(outputs.diagnostics, [
        "testing_unit_id", "role", "event_mass", "status", "event_count",
        "snv_event_count", "indel_event_count",
    ], diagnostics)
    if outputs.partition_model_lookups:
        if len(outputs.partition_model_lookups) != len(params.partition_lookup_keys):
            raise ValueError("Partition/model bucket lookup outputs do not match keys")
        if inputs.partition_statuses is None:
            raise ValueError("Partition/model lookup reduction requires partition statuses")
        pairs = (
            pl.scan_parquet(inputs.partition_statuses)
            .filter(
                (pl.col("status") == "eligible")
                & pl.col("testing_unit_id").is_in(units)
            )
            .select("testing_unit_id", "partition_id", "CHROM")
            .collect(engine="streaming")
            .to_dicts()
        )
        model_rows, model_diagnostics = partition_model_lookup_rows(
            statistics, pairs, hemizygote="hem_chrX" in params.roles, bin_count=params.bin_count
        )
        for path, key in zip(
            outputs.partition_model_lookups, params.partition_lookup_keys, strict=True
        ):
            partition_id, model_id = key.split("\t", 1)
            write_tsv_atomic(path, header, model_rows[(partition_id, model_id)])
        if outputs.partition_model_diagnostics is None:
            raise ValueError("Partition/model lookup reduction requires diagnostics output")
        write_tsv_atomic(
            outputs.partition_model_diagnostics,
            [
                "partition_id", "model", "testing_unit_id", "role", "event_mass",
                "status", "event_count", "snv_event_count", "indel_event_count",
            ],
            ([row[column] for column in (
                "partition_id", "model", "testing_unit_id", "role", "event_mass",
                "status", "event_count", "snv_event_count", "indel_event_count",
            )] for row in model_diagnostics),
        )
