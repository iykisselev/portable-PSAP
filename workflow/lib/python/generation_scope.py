"""Build the shared scope files that every generation window needs."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import polars as pl

from generation_bindings import (
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)
from generation_manifests import (
    BUCKET_MANIFEST_SCHEMA_VERSION,
    WINDOW_MANIFEST_SCHEMA_VERSION,
)
from generation_windows import (
    GenomeWindow,
    build_windows,
    indel_candidate_owner_intervals,
    read_fasta_index,
)
from generation_population import validate_population_af_headers
from generation_prep import (
    effective_calibration_mask,
    merge_intervals,
    overlap_sorted,
    read_bed_units,
    read_intervals,
    supported_scope_mask,
)
from generation_probability import unit_bucket_frame
from partition_audit import partition_id_overlap_rows
from partition_scope import build_partition_scope
from generation_resource_policy import check_unindexed_fixture
from genomic_contracts import GENERATION_CONTIGS
from io_contracts import write_json_atomic, write_parquet_atomic


def sparse_bucket_windows(
    genome_windows: list[GenomeWindow],
    callable_unit_intervals: list[tuple[str, str, int, int]],
    bucket_by_unit: dict[str, int],
    window_size: int,
    *,
    candidate_owner_intervals: list[tuple[str, str, int, int]] | None = None,
) -> dict[int, tuple[str, ...]]:
    """Return only window/bucket pairs with callable candidate overlap.

    ``candidate_owner_intervals`` adds the upstream anchor halo required for
    indels.  It is optional while the active workflow remains SNV-only, so the
    existing sparse topology is byte-for-byte unchanged.
    """
    window_by_coordinate = {
        (window.chromosome, window.start // window_size): window
        for window in genome_windows
    }
    window_ids_by_bucket: dict[int, set[str]] = {
        bucket_id: set() for bucket_id in set(bucket_by_unit.values())
    }
    all_intervals = list(callable_unit_intervals)
    if candidate_owner_intervals is not None:
        all_intervals.extend(candidate_owner_intervals)
    for unit_id, chromosome, start, end in all_intervals:
        bucket_id = bucket_by_unit.get(unit_id)
        if bucket_id is None:
            continue
        first_window = start // window_size
        last_window = (end - 1) // window_size
        for window_index in range(first_window, last_window + 1):
            window = window_by_coordinate.get((chromosome, window_index))
            if window is None:
                raise ValueError(
                    f"Candidate-owner interval for {unit_id!r} has no generation window"
                )
            window_ids_by_bucket[bucket_id].add(window.identifier)
    window_ordinal = {window.identifier: window.ordinal for window in genome_windows}
    return {
        bucket_id: tuple(sorted(identifiers, key=window_ordinal.__getitem__))
        for bucket_id, identifiers in window_ids_by_bucket.items()
    }


UnitInterval = tuple[str, str, int, int]


@dataclass(frozen=True)
class ScopeMasks:
    supported_scope: dict[str, list[tuple[int, int]]]
    calibration_mask: dict[str, list[tuple[int, int]]]


@dataclass(frozen=True)
class ScopeTables:
    units: pl.DataFrame
    unit_intervals: pl.DataFrame
    callable_intervals: list[UnitInterval]
    eligible_intervals: list[UnitInterval]


def _validate_scope_inputs(
    inputs: GenerationScopeInputs, params: GenerationScopeParameters
) -> None:
    # Check the fixture exception before reading resources or scheduling windows.
    for path, index in (
        (inputs.population_vcf, inputs.population_vcf_index),
        (inputs.score_snv, inputs.score_snv_index),
        (inputs.mutation_rate, inputs.mutation_rate_index),
        (
            inputs.score_indel if params.include_indels else None,
            inputs.score_indel_index,
        ),
        (inputs.mutation_rate_indel, inputs.mutation_rate_indel_index),
    ):
        if path is not None and index is None:
            check_unindexed_fixture(path)
    validate_population_af_headers(inputs.population_vcf, dict(params.af_tags))


def _scope_masks(inputs: GenerationScopeInputs) -> ScopeMasks:
    units_by_id: dict[str, list[tuple[str, int, int]]] = {}
    for path in inputs.testing_unit_partitions:
        for unit_id, intervals in read_bed_units(path).items():
            if unit_id in units_by_id and units_by_id[unit_id][0][0] != intervals[0][0]:
                raise ValueError(
                    f"Testing unit {unit_id!r} uses conflicting chromosomes across partitions"
                )
            units_by_id.setdefault(unit_id, []).extend(intervals)
    par_intervals = read_intervals(inputs.par_regions)
    excluded_intervals = read_intervals(inputs.exclusion) if inputs.exclusion else {}
    for exclusion_path in inputs.exclusions:
        for chromosome, intervals in read_intervals(exclusion_path).items():
            excluded_intervals.setdefault(chromosome, []).extend(intervals)

    fasta_index = read_fasta_index(inputs.fasta_fai)
    supported_chromosomes = set(GENERATION_CONTIGS)
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
    return ScopeMasks(supported_scope, calibration_mask)


def _partition_scope_tables(scope, masks: ScopeMasks) -> ScopeTables:
    """Derive shared execution work from eligible resolved ownership only.

    Biological partitions do not receive their own execution buckets.  The
    union of the eligible, calibration-masked owned intervals defines one
    deterministic unit bucket and the sparse window relations for that unit.
    """
    resolved_by_pair = {
        (unit_id, partition_id): intervals
        for partition_id, units in scope.resolved.items()
        for unit_id, intervals in units.items()
    }
    callable_by_unit: dict[str, dict[str, list[tuple[int, int]]]] = {}
    for row in scope.rows:
        if row["status"] != "eligible":
            continue
        unit_id = str(row["testing_unit_id"])
        chromosome = str(row["CHROM"])
        owned = resolved_by_pair[(unit_id, str(row["partition_id"]))]
        retained = overlap_sorted(
            [(start, end) for _, start, end in owned],
            masks.calibration_mask.get(chromosome, []),
        )
        if retained:
            callable_by_unit.setdefault(unit_id, {}).setdefault(chromosome, []).extend(
                retained
            )

    unit_rows: list[tuple[str, str]] = []
    interval_rows: list[tuple[str, str, int, int]] = []
    callable_intervals: list[UnitInterval] = []
    for unit_id, by_chromosome in sorted(callable_by_unit.items()):
        for chromosome, intervals in sorted(by_chromosome.items()):
            for start, end in merge_intervals(intervals):
                interval_rows.append((unit_id, chromosome, start, end))
                callable_intervals.append((unit_id, chromosome, start, end))
                unit_rows.append((unit_id, "eligible"))
    units = pl.DataFrame(
        sorted(set(unit_rows)),
        schema=["testing_unit_id", "status"],
        orient="row",
    )
    unit_intervals = pl.DataFrame(
        interval_rows,
        schema=["testing_unit_id", "CHROM", "START", "END"],
        orient="row",
    ).sort("CHROM", "START", "END", "testing_unit_id")
    return ScopeTables(units, unit_intervals, callable_intervals, callable_intervals)


def _window_manifest(
    genome_windows: list[GenomeWindow], params: GenerationScopeParameters
) -> dict:
    return {
        "schema_version": WINDOW_MANIFEST_SCHEMA_VERSION,
        "window_size_bp": int(params.window_size),
        "windows": [
            {
                "window_id": window.identifier,
                "chromosome": window.chromosome,
                "start": window.start,
                "end": window.end,
                "ordinal": window.ordinal,
            }
            for window in genome_windows
        ],
    }


def _bucket_manifest(
    scope: ScopeTables,
    genome_windows: list[GenomeWindow],
    indel_owner_intervals: list[UnitInterval],
    params: GenerationScopeParameters,
) -> dict:
    unit_retention_frame = scope.units
    unit_intervals_frame = scope.unit_intervals
    callable_unit_intervals = scope.callable_intervals
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
    window_ids_by_bucket = sparse_bucket_windows(
        genome_windows,
        callable_unit_intervals,
        bucket_by_unit,
        int(params.window_size),
        candidate_owner_intervals=indel_owner_intervals
        if params.include_indels
        else None,
    )
    return {
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
                "window_ids": list(window_ids_by_bucket[ordinal]),
            }
            for ordinal, unit_ids in sorted(unit_ids_by_bucket.items())
        ],
    }


def _write_scope_mask(
    masks: ScopeMasks,
    outputs: GenerationScopeOutputs,
    params: GenerationScopeParameters,
) -> None:
    calibration_mask = masks.calibration_mask
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
        outputs.effective_mask,
        params.compression,
        row_group_size=16_384,
    )


def _write_partition_scope_artifacts(
    inputs: GenerationScopeInputs,
    outputs: GenerationScopeOutputs,
    masks: ScopeMasks,
    params: GenerationScopeParameters,
) -> object:
    """Materialize ordered unit/partition ownership and eligibility statuses."""
    paths = inputs.testing_unit_partitions
    identifiers = params.partition_ids
    if len(paths) != len(identifiers):
        raise ValueError(
            "Partition IDs and testing-unit BED inputs must have equal length"
        )
    partitions = [
        (
            identifier,
            {
                unit_id: intervals
                for unit_id, intervals in read_bed_units(path).items()
                if intervals[0][0] in GENERATION_CONTIGS
            },
        )
        for identifier, path in zip(identifiers, paths, strict=True)
    ]
    scope = build_partition_scope(
        partitions,
        masks.supported_scope,
        masks.calibration_mask,
        float(params.minimum_retained_fraction),
    )
    resolved_rows = [
        (unit, partition, order, chromosome, start, end)
        for order, (partition, units) in enumerate(scope.resolved.items())
        for unit, intervals in units.items()
        for chromosome, start, end in intervals
    ]
    write_parquet_atomic(
        pl.DataFrame(
            resolved_rows,
            schema=[
                "testing_unit_id",
                "partition_id",
                "partition_order",
                "CHROM",
                "START",
                "END",
            ],
            orient="row",
        ).sort("CHROM", "START", "END", "testing_unit_id", "partition_id"),
        outputs.resolved_partition_intervals,
        params.compression,
    )
    write_parquet_atomic(
        pl.DataFrame(list(scope.rows)).sort("testing_unit_id", "partition_id"),
        outputs.partition_statuses,
        params.compression,
    )
    audit_rows = partition_id_overlap_rows(partitions)
    schema = {
        "left_partition_id": pl.String,
        "right_partition_id": pl.String,
        "left_unique_id_count": pl.Int64,
        "right_unique_id_count": pl.Int64,
        "intersection_count": pl.Int64,
        "union_count": pl.Int64,
        "left_overlap_fraction": pl.Float64,
        "right_overlap_fraction": pl.Float64,
        "minimum_overlap_fraction": pl.Float64,
        "empty_set_status": pl.String,
        "low_overlap_warning": pl.Boolean,
    }
    write_parquet_atomic(
        pl.DataFrame(audit_rows, schema=schema, strict=False),
        outputs.testing_unit_id_overlap,
        params.compression,
    )
    return scope


def prepare_generation_scope(
    inputs: GenerationScopeInputs,
    outputs: GenerationScopeOutputs,
    params: GenerationScopeParameters,
) -> None:
    """Prepare masks, retention tables and sparse manifests for window jobs."""
    _validate_scope_inputs(inputs, params)
    masks = _scope_masks(inputs)
    partition_scope = _write_partition_scope_artifacts(inputs, outputs, masks, params)
    scope = _partition_scope_tables(partition_scope, masks)
    eligible_callable_unit_intervals = scope.eligible_intervals
    indel_owner_intervals: list[UnitInterval] = []
    if params.include_indels:
        if inputs.score_indel is None:
            raise ValueError("Indel generation requires generation.score_indel")
        indel_owner_intervals = indel_candidate_owner_intervals(
            eligible_callable_unit_intervals,
            params.indel_max_reference_length,
        )
    genome_windows = build_windows(
        inputs.fasta_fai,
        [
            (chromosome, start, end)
            for _, chromosome, start, end in eligible_callable_unit_intervals
        ],
        int(params.window_size),
        candidate_owner_intervals=[
            (chromosome, start, end)
            for _, chromosome, start, end in indel_owner_intervals
        ],
    )
    _write_scope_mask(masks, outputs, params)
    write_json_atomic(outputs.windows, _window_manifest(genome_windows, params))
    write_json_atomic(
        outputs.buckets,
        _bucket_manifest(scope, genome_windows, indel_owner_intervals, params),
    )
