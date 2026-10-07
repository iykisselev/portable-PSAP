"""Build the shared scope files that every generation partition needs."""

from __future__ import annotations

import warnings

import polars as pl

from generation_bindings import (
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)
from generation_manifests import (
    BUCKET_MANIFEST_SCHEMA_VERSION,
    PARTITION_MANIFEST_SCHEMA_VERSION,
)
from generation_partitions import (
    GenomePartition,
    build_partitions,
    indel_candidate_owner_intervals,
    read_fasta_index,
)
from generation_population import validate_population_af_headers
from generation_prep import (
    retention_rows,
    effective_calibration_mask,
    merge_intervals,
    read_bed_units,
    read_intervals,
    supported_scope_mask,
)
from generation_probability import unit_bucket_frame
from generation_resource_policy import check_unindexed_fixture
from io_contracts import write_json_atomic, write_parquet_atomic


def sparse_bucket_partitions(
    genome_partitions: list[GenomePartition],
    callable_unit_intervals: list[tuple[str, str, int, int]],
    bucket_by_unit: dict[str, int],
    partition_size: int,
    *,
    candidate_owner_intervals: list[tuple[str, str, int, int]] | None = None,
) -> dict[int, tuple[str, ...]]:
    """Return only partition/bucket pairs with callable candidate overlap.

    ``candidate_owner_intervals`` adds the upstream anchor halo required for
    indels.  It is optional while the active workflow remains SNV-only, so the
    existing sparse topology is byte-for-byte unchanged.
    """
    partition_by_window = {
        (partition.chromosome, partition.start // partition_size): partition
        for partition in genome_partitions
    }
    partition_ids_by_bucket: dict[int, set[str]] = {
        bucket_id: set() for bucket_id in set(bucket_by_unit.values())
    }
    all_intervals = list(callable_unit_intervals)
    if candidate_owner_intervals is not None:
        all_intervals.extend(candidate_owner_intervals)
    for unit_id, chromosome, start, end in all_intervals:
        bucket_id = bucket_by_unit.get(unit_id)
        if bucket_id is None:
            continue
        first_window = start // partition_size
        last_window = (end - 1) // partition_size
        for window in range(first_window, last_window + 1):
            partition = partition_by_window.get((chromosome, window))
            if partition is None:
                raise ValueError(
                    f"Candidate-owner interval for {unit_id!r} has no generation partition"
                )
            partition_ids_by_bucket[bucket_id].add(partition.identifier)
    partition_ordinal = {
        partition.identifier: partition.ordinal for partition in genome_partitions
    }
    return {
        bucket_id: tuple(sorted(identifiers, key=partition_ordinal.__getitem__))
        for bucket_id, identifiers in partition_ids_by_bucket.items()
    }


def prepare_generation_scope(
    inputs: GenerationScopeInputs,
    outputs: GenerationScopeOutputs,
    params: GenerationScopeParameters,
) -> None:
    """Prepare the shared scope once so partition jobs only read their slice."""
    # Check the fixture exception before reading resources or scheduling partitions.
    for path, index in (
        (inputs.population_vcf, inputs.population_vcf_index),
        (inputs.score_snv, inputs.score_snv_index),
        (inputs.mutation_rate, inputs.mutation_rate_index),
        (inputs.score_indel if params.include_indels else None, inputs.score_indel_index),
        (inputs.mutation_rate_indel, inputs.mutation_rate_indel_index),
    ):
        if path is not None and index is None:
            check_unindexed_fixture(path)
    validate_population_af_headers(inputs.population_vcf, dict(params.af_tags))
    units_by_id = read_bed_units(inputs.testing_units)
    par_intervals = read_intervals(inputs.par_regions)
    excluded_intervals = read_intervals(inputs.exclusion) if inputs.exclusion else {}
    for exclusion_path in inputs.exclusions.values():
        for chromosome, intervals in read_intervals(exclusion_path).items():
            excluded_intervals.setdefault(chromosome, []).extend(intervals)

    fasta_index = read_fasta_index(inputs.fasta_fai)
    supported_chromosomes = {str(value) for value in range(1, 23)} | {"X"}
    unit_chromosomes = {intervals[0][0] for intervals in units_by_id.values()}
    unsupported_unit_chromosomes = unit_chromosomes - supported_chromosomes
    if unsupported_unit_chromosomes:
        ignored = ", ".join(sorted(unsupported_unit_chromosomes))
        warnings.warn(
            "Ignoring testing-unit contigs outside the generated null-distribution "
            f"scope: {ignored}",
            RuntimeWarning,
            stacklevel=2,
        )
        units_by_id = {
            unit_id: intervals
            for unit_id, intervals in units_by_id.items()
            if intervals[0][0] in supported_chromosomes
        }
    supported_scope = supported_scope_mask(
        {
            chromosome: length
            for chromosome, (_, length) in fasta_index.items()
            if chromosome in supported_chromosomes
        },
        par_intervals,
    )
    calibration_mask = effective_calibration_mask(excluded_intervals, supported_scope)
    supported_scope_starts = {
        chromosome: [start for start, _ in intervals]
        for chromosome, intervals in supported_scope.items()
    }
    calibration_mask_starts = {
        chromosome: [start for start, _ in intervals]
        for chromosome, intervals in calibration_mask.items()
    }
    callable_pieces_by_unit: dict[str, list[tuple[int, int]]] = {}
    unit_retention_rows = retention_rows(
        units_by_id,
        calibration_mask,
        float(params.minimum_retained_fraction),
        denominator_scope=supported_scope,
        retained_region_starts=calibration_mask_starts,
        denominator_region_starts=supported_scope_starts,
        callable_intervals_by_unit=callable_pieces_by_unit,
    )
    unit_intervals = []
    callable_unit_intervals = []
    for unit_id, intervals in sorted(units_by_id.items()):
        chromosome = intervals[0][0]
        merged_intervals = merge_intervals(
            [(start, end) for _, start, end in intervals]
        )
        for start, end in merged_intervals:
            unit_intervals.append((unit_id, chromosome, start, end))
        for start, end in callable_pieces_by_unit[unit_id]:
            callable_unit_intervals.append((unit_id, chromosome, start, end))
    unit_retention_frame = pl.DataFrame(
        unit_retention_rows,
        schema=[
            "testing_unit_id",
            "CHROM",
            "supported_bases",
            "retained_bases",
            "retained_fraction",
            "status",
        ],
        orient="row",
    )
    eligible_unit_ids = set(
        unit_retention_frame.filter(pl.col("status") == "eligible")["testing_unit_id"].to_list()
    )
    eligible_callable_unit_intervals = [
        (unit_id, chromosome, start, end)
        for unit_id, chromosome, start, end in callable_unit_intervals
        if unit_id in eligible_unit_ids
    ]
    indel_owner_intervals: list[tuple[str, str, int, int]] = []
    if params.include_indels:
        if inputs.score_indel is None:
            raise ValueError("Indel generation requires generation.score_indel")
        if params.indel_max_reference_length < 1:
            raise ValueError("Indel maximum reference length must be positive")
        indel_owner_intervals = [
            (unit_id, chromosome, start, end)
            for (unit_id, chromosome, _, _), (chromosome, start, end) in zip(
                eligible_callable_unit_intervals,
                indel_candidate_owner_intervals(
                    [
                        (chromosome, start, end)
                        for _, chromosome, start, end in eligible_callable_unit_intervals
                    ],
                    params.indel_max_reference_length,
                ),
            )
        ]
    genome_partitions = build_partitions(
        inputs.fasta_fai,
        [
            (chromosome, start, end)
            for _, chromosome, start, end in eligible_callable_unit_intervals
        ],
        int(params.partition_size),
        candidate_owner_intervals=[
            (chromosome, start, end)
            for _, chromosome, start, end in indel_owner_intervals
        ],
    )
    write_parquet_atomic(unit_retention_frame, outputs.units, params.compression)
    unit_intervals_frame = pl.DataFrame(
        unit_intervals,
        schema=["testing_unit_id", "CHROM", "START", "END"],
        orient="row",
    ).sort("CHROM", "START", "END", "testing_unit_id")
    write_parquet_atomic(
        unit_intervals_frame, outputs.unit_intervals, params.compression,
        row_group_size=16_384,
    )
    write_parquet_atomic(
        pl.DataFrame(
            callable_unit_intervals,
            schema=["testing_unit_id", "CHROM", "START", "END"],
            orient="row",
        ).sort("CHROM", "START", "END", "testing_unit_id"),
        outputs.masked_unit_intervals,
        params.compression,
        row_group_size=16_384,
    )
    write_parquet_atomic(
        pl.DataFrame(
            [
                (chromosome, start, end, "callable_supported_nonexcluded_nonpar")
                for chromosome, intervals in sorted(calibration_mask.items())
                for start, end in intervals
            ],
            schema=["CHROM", "START", "END", "mask_state"],
            orient="row",
        ).sort("CHROM", "START", "END"),
        outputs.mask,
        params.compression,
        row_group_size=16_384,
    )
    write_json_atomic(
        outputs.partitions,
        {
            "schema_version": PARTITION_MANIFEST_SCHEMA_VERSION,
            "partition_size_bp": int(params.partition_size),
            "partitions": [
                {
                    "partition_id": partition.identifier,
                    "chromosome": partition.chromosome,
                    "start": partition.start,
                    "end": partition.end,
                    "ordinal": partition.ordinal,
                }
                for partition in genome_partitions
            ],
        },
    )
    bucket_assignments = unit_bucket_frame(
        unit_retention_frame, int(params.unit_bucket_size), unit_intervals_frame
    )
    bucket_by_unit = {
        str(unit_id): int(bucket_id)
        for unit_id, bucket_id in bucket_assignments.iter_rows()
    }
    unit_ids_by_bucket: dict[int, list[str]] = {}
    for unit_id, bucket_id in bucket_by_unit.items():
        unit_ids_by_bucket.setdefault(int(bucket_id), []).append(str(unit_id))
    partition_ids_by_bucket = sparse_bucket_partitions(
        genome_partitions,
        callable_unit_intervals,
        bucket_by_unit,
        int(params.partition_size),
        candidate_owner_intervals=indel_owner_intervals if params.include_indels else None,
    )
    write_json_atomic(
        outputs.buckets,
        {
            "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
            "unit_bucket_size": int(params.unit_bucket_size),
            "buckets": [
                {
                    "bucket_id": f"{ordinal:05d}",
                    "ordinal": ordinal,
                    "unit_count": len(unit_ids),
                    # Keep manifest identifiers lexically sorted even though
                    # assignment itself follows genomic order.
                    "unit_ids": sorted(unit_ids),
                    "partition_ids": list(partition_ids_by_bucket[ordinal]),
                }
                for ordinal, unit_ids in sorted(unit_ids_by_bucket.items())
            ],
        },
    )
