"""Partition statistics assembled from canonical allele artifacts."""

from __future__ import annotations

import polars as pl

from generation_bindings import (
    PartitionStatisticsInputs,
    PartitionStatisticsOutputs,
    PartitionStatisticsParameters,
)
from generation_event_assembly import assign_effective_frequencies, role_event_frames
from generation_manifests import bucket_from_manifest
from generation_partitions import partition_from_manifest
from generation_probability import statistics_for_role
from generation_scores import ALLELE_COLUMNS, mutation_rate_rows
from io_contracts import write_parquet_atomic


def selected_partition_rows(path: str, candidate_keys: pl.DataFrame) -> pl.DataFrame:
    """Stream only rows whose alleles are present in a bucket's memberships.

    State and score files describe every candidate allele in a genomic
    partition, whereas a calculation job handles only one unit bucket.  A
    semi join keeps the partition files out of the job's resident working set
    and, unlike filtering after an eager read, bounds it by the selected
    bucket's candidate alleles.
    """
    return (
        pl.scan_parquet(path)
        .join(candidate_keys.lazy(), on=ALLELE_COLUMNS, how="semi")
        .collect(engine="streaming")
    )


def candidate_intervals(
    events: pl.DataFrame, variant_class: str
) -> list[tuple[int, int]]:
    """Return merged-source-compatible one-base intervals for candidate alleles."""
    return [
        (int(position) - 1, int(position))
        for (position,) in (
            events.filter(pl.col("variant_class") == variant_class)
            .select("POS")
            .unique()
            .sort("POS")
            .iter_rows()
        )
    ]


def calculate_partition_statistics(
    inputs: PartitionStatisticsInputs,
    outputs: PartitionStatisticsOutputs,
    params: PartitionStatisticsParameters,
    partition_id: str,
    bucket_id: str,
) -> None:
    """Assign absence-policy frequencies and aggregate one partition/bucket."""
    partition = partition_from_manifest(inputs.partitions, partition_id)
    bucket = bucket_from_manifest(inputs.buckets, bucket_id)
    # Filter the high-cardinality membership table before materializing it.
    # A partition may cover many more units than the one bounded bucket being
    # reduced by this job.
    memberships = (
        pl.scan_parquet(inputs.allele_membership)
        .filter(pl.col("testing_unit_id").is_in(bucket.unit_ids))
        .collect(engine="streaming")
    )
    candidate_keys = memberships.select(ALLELE_COLUMNS).unique()
    scores = selected_partition_rows(inputs.allele_scores, candidate_keys)
    states = selected_partition_rows(inputs.allele_states, candidate_keys)
    events = (
        memberships.with_columns(pl.lit(bucket.ordinal).alias("unit_bucket"))
        .join(scores, on=ALLELE_COLUMNS, how="inner", suffix="_score")
        .join(states, on=ALLELE_COLUMNS, how="inner", suffix="_state")
        .filter(pl.col("allele_state") != "rejected")
    )
    class_mismatch = events.filter(
        (pl.col("variant_class") != pl.col("variant_class_score"))
        | (pl.col("variant_class") != pl.col("variant_class_state"))
    )
    if class_mismatch.height:
        raise ValueError(
            "Allele membership, score, and population state disagree on variant class "
            f"in {partition.identifier}"
        )
    events = events.drop("variant_class_score", "variant_class_state")
    if params.unobserved_policy == "mutation_rate":
        snv_intervals = candidate_intervals(events, "snv")
        rates = mutation_rate_rows(
            inputs.mutation_rate,
            inputs.mutation_rate_index,
            partition.chromosome,
            partition.start,
            partition.end,
            intervals=snv_intervals,
        )
        events = events.join(rates, on=ALLELE_COLUMNS, how="left")
        if inputs.mutation_rate_indel is not None:
            indel_intervals = candidate_intervals(events, "indel")
            indel_rates = mutation_rate_rows(
                inputs.mutation_rate_indel,
                inputs.mutation_rate_indel_index,
                partition.chromosome,
                partition.start,
                partition.end,
                intervals=indel_intervals,
            ).rename({"mu": "mu_indel"})
            events = events.join(indel_rates, on=ALLELE_COLUMNS, how="left")
        else:
            events = events.with_columns(
                pl.lit(None, dtype=pl.Float64).alias("mu_indel")
            )
    events = assign_effective_frequencies(
        events,
        unobserved_policy=params.unobserved_policy,
        x_af_policy=params.x_af_policy,
        x_mutation_rate_policy=params.x_mutation_rate_policy,
        roles=params.roles,
    )
    role_events = role_event_frames(events, params.roles)
    statistics = pl.concat(
        [
            statistics_for_role(
                role_events[role], role, params.score_max, params.bin_width
            )
            for role in params.roles
        ],
        how="vertical",
    )
    write_parquet_atomic(
        statistics.drop("unit_bucket"), outputs.statistics, params.compression
    )
