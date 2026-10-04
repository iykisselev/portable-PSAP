"""Build the shared scope files that every generation partition needs."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import polars as pl

from generation_manifests import (
    BUCKET_MANIFEST_SCHEMA_VERSION,
    PARTITION_MANIFEST_SCHEMA_VERSION,
)
from generation_partitions import build_partitions, read_fasta_index
from generation_prep import (
    retention_rows,
    effective_calibration_mask,
    intersect_intervals,
    merge_intervals,
    read_bed_units,
    read_intervals,
    supported_scope_mask,
)
from generation_probability import unit_bucket_frame
from generation_resource_policy import check_unindexed_fixture
from io_contracts import resource_inventory, write_json_atomic, write_parquet_atomic


def prepare_generation_scope(
    inputs: SimpleNamespace, outputs: SimpleNamespace, params: SimpleNamespace
) -> None:
    """Prepare the shared scope once so partition jobs only read their slice."""
    # Check the fixture exception before hashing resources or scheduling partitions.
    for resource_name in ("population_vcf", "score_snv", "mutation_rate"):
        if not getattr(inputs, f"{resource_name}_index", None):
            paths = getattr(inputs, resource_name)
            if isinstance(paths, (str, Path)):
                paths = [paths]
            for path in paths:
                check_unindexed_fixture(path)
    units_by_id = read_bed_units(inputs.testing_units)
    par_intervals = read_intervals(inputs.par_regions)
    excluded_intervals = read_intervals(inputs.exclusion) if inputs.exclusion else {}
    for exclusion_path in inputs.exclusions:
        for chromosome, intervals in read_intervals(exclusion_path).items():
            excluded_intervals.setdefault(chromosome, []).extend(intervals)

    fasta_index = read_fasta_index(inputs.fasta_fai)
    supported_chromosomes = {str(value) for value in range(1, 23)} | {"X"}
    unit_chromosomes = {intervals[0][0] for intervals in units_by_id.values()}
    unsupported_unit_chromosomes = unit_chromosomes - supported_chromosomes
    if unsupported_unit_chromosomes:
        raise ValueError(
            "Testing-unit contigs outside the generated null-distribution scope: "
            + ", ".join(sorted(unsupported_unit_chromosomes))
        )
    supported_scope = supported_scope_mask(
        {
            chromosome: length
            for chromosome, (_, length) in fasta_index.items()
            if chromosome in supported_chromosomes
        },
        par_intervals,
    )
    calibration_mask = effective_calibration_mask(excluded_intervals, supported_scope)
    unit_retention_rows = retention_rows(
        units_by_id,
        calibration_mask,
        float(params.minimum_retained_fraction),
        denominator_scope=supported_scope,
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
        for start, end in intersect_intervals(
            merged_intervals, calibration_mask.get(chromosome, [])
        ):
            callable_unit_intervals.append((unit_id, chromosome, start, end))
    genome_partitions = build_partitions(
        inputs.fasta_fai,
        {chromosome for _, chromosome, _, _ in unit_intervals},
        int(params.partition_size),
    )

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
    write_parquet_atomic(unit_retention_frame, outputs.units, params.compression)
    write_parquet_atomic(
        pl.DataFrame(
            unit_intervals,
            schema=["testing_unit_id", "CHROM", "START", "END"],
            orient="row",
        ),
        outputs.unit_intervals,
        params.compression,
    )
    write_parquet_atomic(
        pl.DataFrame(
            callable_unit_intervals,
            schema=["testing_unit_id", "CHROM", "START", "END"],
            orient="row",
        ),
        outputs.masked_unit_intervals,
        params.compression,
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
        ),
        outputs.mask,
        params.compression,
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
    unit_ids_by_bucket: dict[int, list[str]] = {}
    for unit_id, bucket_id in unit_bucket_frame(
        unit_retention_frame, int(params.unit_bucket_size)
    ).iter_rows():
        unit_ids_by_bucket.setdefault(int(bucket_id), []).append(str(unit_id))
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
                    "unit_ids": unit_ids,
                }
                for ordinal, unit_ids in sorted(unit_ids_by_bucket.items())
            ],
        },
    )
    input_resource_paths = {
        "population_vcf": inputs.population_vcf,
        "score_snv": inputs.score_snv,
        "testing_units": inputs.testing_units,
        "par_regions": inputs.par_regions,
        "fasta": inputs.fasta,
    }
    if inputs.exclusion:
        input_resource_paths["exclusion"] = inputs.exclusion
    input_resource_paths.update(
        {
            f"interval_exclusion_{index:03d}": path
            for index, path in enumerate(inputs.exclusions)
        }
    )
    if inputs.mutation_rate:
        input_resource_paths["mutation_rate"] = inputs.mutation_rate[0]
    write_json_atomic(
        outputs.resource_identities,
        {
            "schema_version": PARTITION_MANIFEST_SCHEMA_VERSION,
            "resources": resource_inventory(input_resource_paths),
        },
    )
