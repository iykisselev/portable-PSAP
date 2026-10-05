"""Partition-local join of reference, score, and population observations."""

from __future__ import annotations

import polars as pl

from generation_bindings import (
    GenerationPartitionInputs,
    GenerationPartitionOutputs,
    GenerationPartitionParameters,
)
from generation_partitions import partition_from_manifest
from generation_population import POPULATION_SCHEMA, population_rows
from generation_prep import intersect_intervals
from generation_reference import ReferenceReader, membership_rows
from generation_scores import ALLELE_COLUMNS, SCORE_SCHEMA, score_rows
from io_contracts import write_parquet_atomic


def prepare_generation_partition(
    inputs: GenerationPartitionInputs,
    outputs: GenerationPartitionOutputs,
    params: GenerationPartitionParameters,
    partition_id: str,
) -> None:
    """Materialize canonical states, memberships, and scores for one partition."""
    partition = partition_from_manifest(inputs.partitions, partition_id)
    units = pl.read_parquet(inputs.units)
    eligible_units = units.filter(pl.col("status") == "eligible").select("testing_unit_id").sort("testing_unit_id")
    intervals = (
        pl.read_parquet(inputs.unit_intervals)
        .filter((pl.col("CHROM") == partition.chromosome) & (pl.col("START") < partition.end) & (pl.col("END") > partition.start))
        .join(eligible_units, on="testing_unit_id", how="inner")
        .select("testing_unit_id", "START", "END")
        .rows()
    )
    mask = (
        pl.read_parquet(inputs.mask)
        .filter((pl.col("CHROM") == partition.chromosome) & (pl.col("START") < partition.end) & (pl.col("END") > partition.start))
        .select(pl.max_horizontal("START", pl.lit(partition.start)).alias("START"), pl.min_horizontal("END", pl.lit(partition.end)).alias("END"))
        .rows()
    )
    reference = ReferenceReader(inputs.fasta, inputs.fasta_fai)
    try:
        memberships = membership_rows(partition, mask, intervals, reference)
    finally:
        reference.close()
    candidates = memberships.select(ALLELE_COLUMNS).unique(maintain_order=True)
    if candidates.is_empty():
        # Keep the same artifact contract without opening either large resource.
        write_parquet_atomic(pl.DataFrame(schema=POPULATION_SCHEMA), outputs.allele_states, params.compression)
        write_parquet_atomic(memberships, outputs.allele_membership, params.compression)
        write_parquet_atomic(pl.DataFrame(schema=SCORE_SCHEMA), outputs.allele_scores, params.compression)
        return
    candidate_intervals = intersect_intervals([(start, end) for _, start, end in intervals], mask)
    scores = score_rows(
        inputs.score_snv, inputs.score_snv_index, partition.chromosome, partition.start, partition.end,
        intervals=candidate_intervals,
    )
    matched_score_keys = candidates.join(scores.select(ALLELE_COLUMNS), on=ALLELE_COLUMNS, how="inner")
    if matched_score_keys.select(ALLELE_COLUMNS).is_duplicated().any():
        raise ValueError(f"Score resource has duplicate calibration alleles in {partition.identifier}")
    scored = candidates.join(scores, on=ALLELE_COLUMNS, how="inner")
    if scored.height != candidates.height:
        missing = candidates.join(scores.select(ALLELE_COLUMNS), on=ALLELE_COLUMNS, how="anti")
        examples = ", ".join(f"{row[0]}:{row[1]}:{row[2]}:{row[3]}" for row in missing.head(10).iter_rows())
        raise ValueError(f"Score resource is missing {missing.height} calibration alleles in {partition.identifier}: {examples}")
    panel = population_rows(
        inputs.population_vcf, inputs.population_vcf_index, partition.chromosome, partition.start,
        partition.end, dict(params.af_tags), params.x_af_policy, params.hemizygote,
    )
    if panel.select(ALLELE_COLUMNS).is_duplicated().any():
        raise ValueError(f"Population VCF has duplicate allele states in {partition.identifier}")
    states = (
        candidates.join(panel, on=ALLELE_COLUMNS, how="left")
        .with_columns(
            pl.when(pl.col("allele_state").is_null()).then(pl.lit("absent")).otherwise(pl.col("allele_state")).alias("allele_state"),
            pl.when(pl.col("state_reason").is_null()).then(pl.lit("not_present_in_population_panel")).otherwise(pl.col("state_reason")).alias("state_reason"),
        )
        .select(
            *ALLELE_COLUMNS, "allele_state", "af", "af_overall", "af_female",
            "af_male", "state_reason", "af_zero_imputed_tags",
        )
    )
    invalid = states.filter(pl.col("allele_state") == "invalid")
    if invalid.height:
        examples = "; ".join(
            f"{chrom}:{position}:{reference}:{alternate}: {reason}"
            for chrom, position, reference, alternate, reason in invalid.head(10).select(*ALLELE_COLUMNS, "state_reason").iter_rows()
        )
        raise ValueError(
            f"Population VCF has invalid selected AF for {invalid.height} calibration allele(s) in {partition.identifier}: {examples}"
        )
    write_parquet_atomic(states, outputs.allele_states, params.compression)
    write_parquet_atomic(memberships, outputs.allele_membership, params.compression)
    write_parquet_atomic(scored, outputs.allele_scores, params.compression)
