"""Partition statistics assembled from canonical allele artifacts."""

from __future__ import annotations

import polars as pl

from generation_bindings import (
    PartitionStatisticsInputs,
    PartitionStatisticsOutputs,
    PartitionStatisticsParameters,
)
from generation_manifests import bucket_from_manifest
from generation_partitions import partition_from_manifest
from generation_probability import statistics_for_role
from generation_scores import ALLELE_COLUMNS, mutation_rate_rows
from io_contracts import write_parquet_atomic


def calculate_partition_statistics(
    inputs: PartitionStatisticsInputs,
    outputs: PartitionStatisticsOutputs,
    params: PartitionStatisticsParameters,
    partition_id: str,
    bucket_id: str,
) -> None:
    """Assign absence-policy frequencies and aggregate one partition/bucket."""
    partition = partition_from_manifest(inputs.partitions, partition_id)
    states = pl.read_parquet(inputs.allele_states)
    memberships = pl.read_parquet(inputs.allele_membership)
    scores = pl.read_parquet(inputs.allele_scores)
    bucket = bucket_from_manifest(inputs.buckets, bucket_id)
    selected_units = pl.DataFrame({"testing_unit_id": list(bucket.unit_ids)})
    events = (
        memberships.join(selected_units, on="testing_unit_id", how="inner")
        .with_columns(pl.lit(bucket.ordinal).alias("unit_bucket"))
        .join(scores, on=ALLELE_COLUMNS, how="inner")
        .join(states, on=ALLELE_COLUMNS, how="inner")
        .filter(pl.col("allele_state") != "rejected")
    )
    if params.x_af_policy not in {"sex_specific", "overall_proxy"}:
        raise ValueError(f"Unsupported X chromosome AF policy {params.x_af_policy!r}")
    female_af = pl.when(pl.col("CHROM") == "X").then(
        pl.col("af_female") if params.x_af_policy == "sex_specific" else pl.col("af_overall")
    ).otherwise(pl.col("af_overall"))
    male_af = pl.col("af_male") if params.x_af_policy == "sex_specific" else pl.col("af_overall")
    if params.unobserved_policy == "zero":
        events = events.with_columns(
            pl.when(pl.col("allele_state") == "absent").then(pl.lit(0.0)).otherwise(female_af).alias("effective_af_standard"),
            pl.when(pl.col("allele_state") == "absent").then(pl.lit(0.0)).otherwise(male_af).alias("effective_af_hem"),
        )
    else:
        rates = mutation_rate_rows(inputs.mutation_rate, inputs.mutation_rate_index, partition.chromosome, partition.start, partition.end)
        if params.x_mutation_rate_policy == "reject" and events.filter((pl.col("CHROM") == "X") & (pl.col("allele_state") == "absent")).height:
            raise ValueError(
                "mutation_rate fallback encountered absent chrX alleles, but generation_settings.x_chromosome.mutation_rate_policy is reject"
            )
        events = events.join(rates, on=ALLELE_COLUMNS, how="left")
        missing_rates = events.filter((pl.col("allele_state") == "absent") & pl.col("mu").is_null())
        if missing_rates.height:
            raise ValueError(f"Mutation-rate resource is missing {missing_rates.height} absent alleles in {partition.identifier}")
        events = events.with_columns(
            pl.when(pl.col("allele_state") == "absent").then(pl.col("mu")).otherwise(female_af).alias("effective_af_standard"),
            pl.when(pl.col("allele_state") == "absent").then(pl.col("mu")).otherwise(male_af).alias("effective_af_hem"),
        )
    role_events = {
        "het": events.with_columns(pl.col("effective_af_standard").alias("effective_af")),
        "hom": events.with_columns(pl.col("effective_af_standard").alias("effective_af")),
        "hem_chrX": events.filter(pl.col("CHROM") == "X").with_columns(pl.col("effective_af_hem").alias("effective_af")),
    }
    statistics = pl.concat(
        [statistics_for_role(role_events[role], role, params.score_max, params.bin_width) for role in params.roles],
        how="vertical",
    )
    write_parquet_atomic(statistics.drop("unit_bucket"), outputs.statistics, params.compression)
