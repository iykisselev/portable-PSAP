"""Render partition/model lookup rows from partition-keyed role statistics."""

from __future__ import annotations

from collections import defaultdict

import polars as pl

from generation_probability import vectorized_lookup_rows
from model_catalog import DiseaseModel, enabled_models


def _model_units(rows: list[dict[str, object]], model: DiseaseModel) -> list[str]:
    chromosome = "X" if model.scope == "non_PAR_chrX" else None
    return [
        str(row["testing_unit_id"])
        for row in rows
        if (str(row["CHROM"]) == chromosome if chromosome else str(row["CHROM"]) != "X")
    ]


def partition_model_lookup_rows(
    statistics: pl.DataFrame,
    eligible_pairs: list[dict[str, object]],
    *,
    hemizygote: bool,
    bin_count: int,
) -> tuple[dict[tuple[str, str], list[tuple[object, ...]]], list[dict[str, object]]]:
    """Return separate rows for every configured partition and enabled model.

    ``eligible_pairs`` supplies the model scope even when a partition/model has
    no events, so callers can publish a valid header-only physical table.
    """
    required = {"testing_unit_id", "partition_id", "role", "bin"}
    missing = required - set(statistics.columns)
    if missing:
        raise ValueError("Partition statistics are missing columns: " + ", ".join(sorted(missing)))
    by_partition: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in eligible_pairs:
        by_partition[str(row["partition_id"])].append(row)
    result: dict[tuple[str, str], list[tuple[object, ...]]] = {}
    diagnostics: list[dict[str, object]] = []
    for partition_id, pairs in by_partition.items():
        partition_statistics = statistics.filter(pl.col("partition_id") == partition_id)
        for model in enabled_models(hemizygote):
            units = _model_units(pairs, model)
            role_statistics = (
                partition_statistics.filter(
                    (pl.col("role") == model.role)
                    & pl.col("testing_unit_id").is_in(units)
                )
                .drop("partition_id")
            )
            lookups, model_diagnostics = vectorized_lookup_rows(
                role_statistics, units, (model.role,), bin_count
            )
            result[(partition_id, model.identifier)] = lookups[model.role]
            diagnostics.extend(
                {
                    "partition_id": partition_id,
                    "model": model.identifier,
                    "testing_unit_id": row[0],
                    "role": row[1],
                    "event_mass": row[2],
                    "status": row[3],
                    "event_count": row[4],
                    "snv_event_count": row[5],
                    "indel_event_count": row[6],
                }
                for row in model_diagnostics
            )
    return result, diagnostics
