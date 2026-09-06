"""Construction of data-dependent generation scope artifacts."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import polars as pl

try:
    from generation_manifests import BUCKET_MANIFEST_SCHEMA_VERSION, PARTITION_MANIFEST_SCHEMA_VERSION
    from generation_partitions import build_partitions, read_fasta_index
    from generation_prep import coverage_rows, effective_calibration_mask, intersect_intervals, merge_intervals, read_bed_units, read_intervals, supported_scope_mask
    from generation_probability import unit_bucket_frame
    from io_contracts import resource_inventory, write_json_atomic, write_parquet_atomic
except ModuleNotFoundError:
    from .generation_manifests import BUCKET_MANIFEST_SCHEMA_VERSION, PARTITION_MANIFEST_SCHEMA_VERSION
    from .generation_partitions import build_partitions, read_fasta_index
    from .generation_prep import coverage_rows, effective_calibration_mask, intersect_intervals, merge_intervals, read_bed_units, read_intervals, supported_scope_mask
    from .generation_probability import unit_bucket_frame
    from .io_contracts import resource_inventory, write_json_atomic, write_parquet_atomic


def prepare_generation_scope(
    inputs: SimpleNamespace, outputs: SimpleNamespace, params: SimpleNamespace
) -> None:
    """Normalize static scope once; later jobs consume only bounded artifacts."""
    units = read_bed_units(inputs.testing_units)
    coverage = read_intervals(inputs.coverage)
    exclusions = read_intervals(inputs.par_regions)
    for exclusion_path in inputs.exclusions:
        for chromosome, intervals in read_intervals(exclusion_path).items():
            exclusions[chromosome].extend(intervals)

    fasta_index = read_fasta_index(inputs.fasta_fai)
    supported_chromosomes = {str(value) for value in range(1, 23)} | {"X"}
    unit_chromosomes = {intervals[0][0] for intervals in units.values()}
    unsupported_unit_chromosomes = unit_chromosomes - supported_chromosomes
    if unsupported_unit_chromosomes:
        raise ValueError(
            "Testing-unit contigs outside the generated null-distribution scope: "
            + ", ".join(sorted(unsupported_unit_chromosomes))
        )
    supported_scope = supported_scope_mask(
        {chromosome: length for chromosome, (_, length) in fasta_index.items() if chromosome in supported_chromosomes},
        exclusions,
    )
    mask = effective_calibration_mask(coverage, supported_scope)
    unit_rows = coverage_rows(units, mask, float(params.coverage_minimum), denominator_scope=supported_scope)
    normalized_intervals = []
    masked_intervals = []
    for unit, intervals in sorted(units.items()):
        chromosome = intervals[0][0]
        merged = merge_intervals([(start, end) for _, start, end in intervals])
        for start, end in merged:
            normalized_intervals.append((unit, chromosome, start, end))
        for start, end in intersect_intervals(merged, mask.get(chromosome, [])):
            masked_intervals.append((unit, chromosome, start, end))
    partitions = build_partitions(inputs.fasta_fai, {chromosome for _, chromosome, _, _ in normalized_intervals}, int(params.partition_size))

    unit_frame = pl.DataFrame(
        unit_rows,
        schema=["testing_unit_id", "CHROM", "union_bases", "covered_bases", "coverage_fraction", "status"],
        orient="row",
    )
    write_parquet_atomic(unit_frame, outputs.units, params.compression)
    write_parquet_atomic(
        pl.DataFrame(normalized_intervals, schema=["testing_unit_id", "CHROM", "START", "END"], orient="row"),
        outputs.unit_intervals,
        params.compression,
    )
    write_parquet_atomic(
        pl.DataFrame(masked_intervals, schema=["testing_unit_id", "CHROM", "START", "END"], orient="row"),
        outputs.masked_unit_intervals,
        params.compression,
    )
    write_parquet_atomic(
        pl.DataFrame(
            [(chromosome, start, end, "callable_supported_nonexcluded_nonpar") for chromosome, intervals in sorted(mask.items()) for start, end in intervals],
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
                {"partition_id": partition.identifier, "chromosome": partition.chromosome, "start": partition.start, "end": partition.end, "ordinal": partition.ordinal}
                for partition in partitions
            ],
        },
    )
    bucket_units: dict[int, list[str]] = {}
    for unit, bucket in unit_bucket_frame(unit_frame, int(params.unit_bucket_size)).iter_rows():
        bucket_units.setdefault(int(bucket), []).append(str(unit))
    write_json_atomic(
        outputs.buckets,
        {
            "schema_version": BUCKET_MANIFEST_SCHEMA_VERSION,
            "unit_bucket_size": int(params.unit_bucket_size),
            "buckets": [
                {"bucket_id": f"{ordinal:05d}", "ordinal": ordinal, "unit_count": len(unit_ids), "unit_ids": unit_ids}
                for ordinal, unit_ids in sorted(bucket_units.items())
            ],
        },
    )
    resource_paths = {
        "population_vcf": inputs.population_vcf,
        "score_snv": inputs.score_snv,
        "coverage": inputs.coverage,
        "testing_units": inputs.testing_units,
        "par_regions": inputs.par_regions,
        "fasta": inputs.fasta,
    }
    resource_paths.update({f"interval_exclusion_{index:03d}": path for index, path in enumerate(inputs.exclusions)})
    if inputs.mutation_rate:
        resource_paths["mutation_rate"] = inputs.mutation_rate[0]
    write_json_atomic(
        outputs.resource_identities,
        {"schema_version": PARTITION_MANIFEST_SCHEMA_VERSION, "resources": resource_inventory(resource_paths)},
    )
