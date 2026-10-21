"""Streaming assembly of generated null-distribution bundle components."""

from __future__ import annotations

import csv
from collections import defaultdict
import json
from pathlib import Path

import polars as pl

from generation_bindings import (
    BundleFinalizationInputs,
    BundleFinalizationOutputs,
    BundleFinalizationParameters,
)
from bundle_contract import (
    CANONICAL_CONTIG_CONVENTION,
    CURRENT_GENERATED_BUNDLE_CONTRACT,
    CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION,
    FIXED_WIDTH_BINNING_ALGORITHM,
    VALIDATION_ARTIFACT_FILENAMES,
    model_roles,
)
from io_contracts import GzipTsvAtomicWriter, write_gzip_tsv_atomic, write_json_atomic
from model_catalog import enabled_models, models_for_chromosome
from partitioned_bundle_contract import (
    PARTITIONED_NULL_BUNDLE_CONTRACT,
    PARTITIONED_NULL_BUNDLE_SCHEMA_VERSION,
)


def lookup_rows(paths, allowed_units):
    """Copy direct bucket lookup rows without constructing a long lookup table."""
    for path in paths:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader, None)
            for row in reader:
                if row and row[0] in allowed_units:
                    yield row


def parquet_rows(paths, columns: list[str]):
    """Yield selected Parquet columns in record batches rather than whole tables."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        for frame in parquet_frames(paths, columns):
            yield from frame.iter_rows()
        return
    for path in paths:
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(columns=columns, batch_size=8192):
            values = batch.to_pydict()
            yield from zip(*(values[column] for column in columns), strict=True)


def parquet_frames(paths, columns: list[str]):
    """Yield small selected Parquet frames without materializing a whole shard."""
    for path in paths:
        yield from (
            pl.scan_parquet(path)
            .select(columns)
            .collect_batches(chunk_size=8192, engine="streaming")
        )


def candidate_audit_rows(paths):
    """Expose structural indel candidate rejections for the public audit."""
    for chromosome, position, reference, alternate, reason in parquet_rows(
        paths, ["CHROM", "POS", "REF", "ALT", "rejection_reason"]
    ):
        yield chromosome, position, reference, alternate, reason


def _empty_variant_class_counts() -> dict[str, dict[str, int]]:
    return {
        variant_class: {
            name: 0 for name in ("retained", "rejected", "observed", "absent")
        }
        for variant_class in ("snv", "indel")
    }


def _update_variant_class_counts(
    counts: dict[str, dict[str, int]],
    variant_class: str,
    state: str,
    amount: int = 1,
) -> None:
    """Add one grouped population-state count to the v3 metadata counts."""
    if variant_class not in counts:
        raise ValueError(
            f"Unsupported variant class in allele states: {variant_class!r}"
        )
    if state in {"observed", "absent"}:
        counts[variant_class]["retained"] += amount
        counts[variant_class][state] += amount
    elif state == "rejected":
        counts[variant_class]["rejected"] += amount
    else:
        raise ValueError(f"Unsupported allele state in bundle finalization: {state!r}")


def _variant_class_counts_from_frame(
    count_frame: pl.DataFrame,
) -> dict[str, dict[str, int]]:
    """Build public v3 counts from a compact grouped state-count frame."""
    counts = _empty_variant_class_counts()
    for variant_class, state, amount in count_frame.iter_rows():
        _update_variant_class_counts(counts, variant_class, state, amount)
    return counts


def _merge_variant_class_counts(
    target: dict[str, dict[str, int]], source: dict[str, dict[str, int]]
) -> None:
    """Add compact per-shard counts without retaining a candidate-state frame."""
    for variant_class, source_counts in source.items():
        for state, amount in source_counts.items():
            target[variant_class][state] += amount


def variant_class_counts(state_paths, audit_paths) -> dict[str, dict[str, int]]:
    """Count retained, population-rejected, and structural candidates by class."""
    counts = _empty_variant_class_counts()
    for state_path in state_paths:
        count_frame = (
            pl.scan_parquet(state_path)
            .group_by("variant_class", "allele_state")
            .len()
            .collect(engine="streaming")
        )
        _merge_variant_class_counts(
            counts, _variant_class_counts_from_frame(count_frame)
        )
    for _ in candidate_audit_rows(audit_paths):
        counts["indel"]["rejected"] += 1
    return counts


def write_finalization_allele_artifacts(
    state_paths,
    audit_paths,
    outputs: BundleFinalizationOutputs,
    *,
    include_variant_class_counts: bool,
) -> tuple[
    int,
    dict[str, object],
    dict[str, object],
    int,
    dict[str, object],
    dict[str, dict[str, int]] | None,
]:
    """Publish sparse state artifacts in one bounded pass over every shard.

    Each state batch is capped at 8,192 rows. The three gzip outputs stay open
    concurrently, so state rows are read once and audit rows are not retained
    merely to copy them into two public artifacts. Excluded allele ordering is
    still population-state rejections followed by structural indel rejections.
    """
    counts = _empty_variant_class_counts() if include_variant_class_counts else None
    state_columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "allele_state",
        "state_reason",
        "af_zero_imputed_tags",
    ]
    if counts is not None:
        state_columns.append("variant_class")
    with (
        GzipTsvAtomicWriter(
            outputs.excluded_alleles,
            ["#CHROM", "POS", "REF", "ALT", "exclusion_reason"],
        ) as excluded,
        GzipTsvAtomicWriter(
            outputs.candidate_audit,
            ["#CHROM", "POS", "REF", "ALT", "rejection_reason"],
        ) as audit,
        GzipTsvAtomicWriter(
            outputs.population_af_zero_imputations,
            ["#CHROM", "POS", "REF", "ALT", "imputed_info_tags"],
        ) as zero_imputed,
    ):
        for frame in parquet_frames(state_paths, state_columns):
            if counts is not None:
                count_frame = frame.group_by("variant_class", "allele_state").len()
                _merge_variant_class_counts(
                    counts, _variant_class_counts_from_frame(count_frame)
                )
            for row in (
                frame.filter(pl.col("allele_state") == "rejected")
                .select("CHROM", "POS", "REF", "ALT", "state_reason")
                .iter_rows()
            ):
                excluded.write_row(row)
            for row in (
                frame.filter(
                    pl.col("af_zero_imputed_tags").is_not_null()
                    & (pl.col("af_zero_imputed_tags") != "")
                )
                .select("CHROM", "POS", "REF", "ALT", "af_zero_imputed_tags")
                .iter_rows()
            ):
                zero_imputed.write_row(row)
        for frame in parquet_frames(
            audit_paths, ["CHROM", "POS", "REF", "ALT", "rejection_reason"]
        ):
            for row in frame.iter_rows():
                excluded.write_row(row)
                audit.write_row(row)
                if counts is not None:
                    counts["indel"]["rejected"] += 1
    assert (
        excluded.inventory is not None
        and audit.inventory is not None
        and zero_imputed.inventory is not None
    )
    return (
        excluded.count,
        excluded.inventory,
        audit.inventory,
        zero_imputed.count,
        zero_imputed.inventory,
        counts,
    )


def finalize_generated_bundle(
    inputs: BundleFinalizationInputs,
    outputs: BundleFinalizationOutputs,
    params: BundleFinalizationParameters,
) -> None:
    """Build a certified public bundle from bounded generation artifacts."""
    combined_indel_bundle = params.variant_classes == ("snv", "indel")
    if params.variant_classes not in {("snv",), ("snv", "indel")}:
        raise ValueError(
            f"Unsupported generated variant classes: {params.variant_classes!r}"
        )
    units = pl.read_parquet(inputs.units).to_dicts()
    eligible, excluded = classify_units(
        units, list(diagnostics_rows(inputs.diagnostics)), hemizygote=params.hemizygote
    )
    eligible_set = set(eligible)
    eligible_x_units = {
        row["testing_unit_id"]
        for row in units
        if row["CHROM"] == "X" and row["testing_unit_id"] in eligible_set
    }
    if params.hemizygote and not eligible_x_units:
        raise ValueError(
            "HEM generation requested, but no chromosome-X testing unit survived the HET/HOM/HEM eligibility contract"
        )
    header = [
        "testing_unit_id",
        *[str(value) for value in range(1, params.bin_count + 1)],
    ]
    _, lookup_het = write_gzip_tsv_atomic(
        outputs.lookup_het, header, lookup_rows(inputs.lookup_het, eligible_set)
    )
    _, lookup_hom = write_gzip_tsv_atomic(
        outputs.lookup_hom, header, lookup_rows(inputs.lookup_hom, eligible_set)
    )
    lookup_hem_chrX = None
    if params.hemizygote:
        hem_count, lookup_hem_chrX = write_gzip_tsv_atomic(
            outputs.lookup_hem_chrX,
            header,
            lookup_rows(inputs.lookup_hem_chrX, eligible_x_units),
        )
        if hem_count != len(eligible_x_units):
            raise ValueError(
                "HEM lookup rows do not exactly cover surviving chromosome-X testing units"
            )
    expected_lookup_count = len(params.partition_ids) * len(enabled_models(params.hemizygote))
    if outputs.partition_lookups and len(outputs.partition_lookups) != expected_lookup_count:
        raise ValueError("Partition/model lookup output count does not match the bundle layout")
    partition_components = {}
    if outputs.partition_lookups:
        if len(params.partition_ids) != 1:
            raise ValueError("Partition-aware finalization is required for N > 1")
        partition_id = params.partition_ids[0]
        role_inputs = {
            "het": inputs.lookup_het,
            "hom": inputs.lookup_hom,
            "hem_chrX": inputs.lookup_hem_chrX,
        }
        model_outputs = iter(outputs.partition_lookups)
        for model in enabled_models(params.hemizygote):
            allowed_units = eligible_x_units if model.scope == "non_PAR_chrX" else (
                eligible_set - eligible_x_units
            )
            _, inventory = write_gzip_tsv_atomic(
                next(model_outputs),
                header,
                lookup_rows(role_inputs[model.role], allowed_units),
            )
            # Keep the descriptor materialized in metadata until the complete
            # partitioned certificate replaces the role-based contract.
            partition_components[f"lookup.{partition_id}.{model.identifier}"] = inventory
    _, excluded_units = write_gzip_tsv_atomic(
        outputs.excluded_units, ["testing_unit_id", "exclusion_reason"], excluded
    )
    _, eligible_bed = write_gzip_tsv_atomic(
        outputs.eligible_bed,
        ["#CHROM", "START", "END", "testing_unit_id"],
        masked_bed_rows(inputs.masked_unit_intervals, eligible_set),
    )
    calibration_count, calibration_regions = write_gzip_tsv_atomic(
        outputs.calibration_regions,
        ["#CHROM", "START", "END", "source"],
        parquet_rows([inputs.effective_mask], ["CHROM", "START", "END", "mask_state"]),
    )
    (
        rejected_count,
        excluded_alleles,
        candidate_audit,
        zero_imputation_count,
        population_af_zero_imputations,
        class_counts,
    ) = write_finalization_allele_artifacts(
        inputs.allele_states,
        inputs.candidate_audit,
        outputs,
        include_variant_class_counts=combined_indel_bundle,
    )
    components = {
        "lookup_het": lookup_het,
        "lookup_hom": lookup_hom,
        "eligible_testing_units": eligible_bed,
        "calibration_regions": calibration_regions,
        "excluded_alleles": excluded_alleles,
        "candidate_audit": candidate_audit,
        "population_af_zero_imputations": population_af_zero_imputations,
        "excluded_units": excluded_units,
    }
    if params.hemizygote:
        assert lookup_hem_chrX is not None
        components["lookup_hem_chrX"] = lookup_hem_chrX
    components.update(partition_components)
    metadata = {
        "schema_version": (
            INDEL_GENERATED_BUNDLE_SCHEMA_VERSION
            if combined_indel_bundle
            else CURRENT_GENERATED_BUNDLE_SCHEMA_VERSION
        ),
        "contract": (
            INDEL_GENERATED_BUNDLE_CONTRACT
            if combined_indel_bundle
            else CURRENT_GENERATED_BUNDLE_CONTRACT
        ),
        "assembly": params.assembly,
        "contig_convention": CANONICAL_CONTIG_CONVENTION,
        "testing_unit_type": params.testing_unit_type,
        "af_info_tags": dict(params.af_tags),
        "binning": {
            "algorithm": FIXED_WIDTH_BINNING_ALGORITHM,
            "score_max": params.score_max,
            "bin_width": params.bin_width,
            "bin_count": params.bin_count,
            "maximum_in_final_bin": True,
        },
        "model_assumptions": {
            "minimum_retained_fraction": params.minimum_retained_fraction,
            "retention_denominator": "supported_non_PAR_unit_bases_before_interval_exclusions",
            "roles": list(model_roles(params.hemizygote)),
            "variant_classes": list(params.variant_classes),
            "unobserved_policy": params.policy,
            "mutation_conversion": params.conversion
            if params.policy == "mutation_rate"
            else None,
            "allele_events_independent": True,
            "par_policy": "excluded_before_event_assignment",
            "hemizygote_role": "hem_chrX" if params.hemizygote else None,
            "x_chromosome": {
                "af_policy": params.x_af_policy,
                "mutation_rate_policy": params.x_mutation_rate_policy,
                "hemizygote_scope": "non_PAR_chrX_only",
            },
        },
        "counts": {
            "eligible_units": len(eligible),
            "excluded_units": len(excluded),
            "excluded_alleles": rejected_count,
            "calibration_intervals": calibration_count,
            "population_af_zero_imputations": zero_imputation_count,
        },
        "resource_identities": dict(params.resource_identities),
        "components": components,
        "validation_artifacts": dict(VALIDATION_ARTIFACT_FILENAMES),
    }
    if combined_indel_bundle:
        metadata["indel_universe"] = {
            "candidate_universe": params.indel_candidate_universe,
            "normalization": "vcf_minimal_left_aligned_v1",
            "reference_span": "zero_based_half_open_ref_span_v1",
            "membership": params.indel_membership,
            "mask_policy": params.indel_mask_policy,
            "max_reference_length": params.indel_max_reference_length,
            "max_alternate_length": params.indel_max_alternate_length,
        }
        assert class_counts is not None
        metadata["variant_class_counts"] = class_counts
    write_json_atomic(
        outputs.metadata,
        metadata,
    )


def _partition_model_bucket_paths(
    paths: tuple[str | Path, ...], keys: tuple[str, ...]
) -> dict[tuple[str, str], tuple[str | Path, ...]]:
    if not keys or len(paths) % len(keys):
        raise ValueError("Partition/model lookup bucket paths do not match declared keys")
    grouped: dict[tuple[str, str], list[str | Path]] = defaultdict(list)
    for offset, path in enumerate(paths):
        partition_id, model_id = keys[offset % len(keys)].split("\t", 1)
        grouped[(partition_id, model_id)].append(path)
    return {key: tuple(value) for key, value in grouped.items()}


def _partition_model_diagnostic_rows(paths):
    for path in paths:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def _partition_statuses_for_publication(inputs, params):
    if inputs.partition_statuses is None:
        raise ValueError("Partitioned bundle finalization requires pair statuses")
    status_rows = pl.read_parquet(inputs.partition_statuses).to_dicts()
    diagnostics: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in _partition_model_diagnostic_rows(inputs.partition_model_diagnostics):
        diagnostics[(row["testing_unit_id"], row["partition_id"], row["model"])].append(row)
    published: list[dict[str, object]] = []
    for row in status_rows:
        failed_roles: list[str] = []
        if row["status"] == "eligible":
            for model in models_for_chromosome(row["CHROM"], params.hemizygote):
                values = diagnostics.get(
                    (row["testing_unit_id"], row["partition_id"], model.identifier)
                )
                if not values:
                    raise ValueError(
                        "Missing partition/model diagnostic for "
                        f"{row['testing_unit_id']}/{row['partition_id']}/{model.identifier}"
                    )
                if len(values) != 1:
                    raise ValueError("Duplicate partition/model diagnostic")
                value = values[0]
                validate_status(value["status"], "role", role=model.role)
                if value["status"] == "numerical_failure":
                    raise ValueError("Numerical distribution failure prevents bundle publication")
                if value["status"] != "ok":
                    failed_roles.append(model.role)
            if failed_roles:
                row["status"] = "zero_event_mass"
        row["failed_roles_json"] = json.dumps(
            [role for role in GENERATED_ROLES if role in set(failed_roles)], separators=(",", ":")
        )
        published.append(row)
    return published


def _write_partition_bed(path, rows):
    return write_gzip_tsv_atomic(path, ["#CHROM", "START", "END", "testing_unit_id"], rows)


def finalize_partitioned_bundle(
    inputs: BundleFinalizationInputs,
    outputs: BundleFinalizationOutputs,
    params: BundleFinalizationParameters,
) -> None:
    """Publish the sole partitioned-null bundle contract from bounded shards."""
    if not params.partition_ids or tuple(dict.fromkeys(params.partition_ids)) != params.partition_ids:
        raise ValueError("Partitioned bundle requires ordered unique partition IDs")
    expected_keys = tuple(
        f"{partition_id}\t{model.identifier}"
        for partition_id in params.partition_ids
        for model in enabled_models(params.hemizygote)
    )
    if params.partition_lookup_keys != expected_keys:
        raise ValueError("Partition/model lookup keys do not match the enabled model catalog")
    if len(outputs.partition_lookups) != len(expected_keys):
        raise ValueError("Partition/model final lookup outputs do not match declared keys")
    if len(outputs.resolved_partitions) != len(params.partition_ids) or len(outputs.eligible_partitions) != len(params.partition_ids):
        raise ValueError("Partition BED output count does not match configured partitions")
    if outputs.partition_statuses is None or outputs.testing_unit_id_overlap is None:
        raise ValueError("Partitioned bundle requires public scope audit outputs")

    status_rows = _partition_statuses_for_publication(inputs, params)
    eligible_pairs = {
        (str(row["testing_unit_id"]), str(row["partition_id"]))
        for row in status_rows if row["status"] == "eligible"
    }
    lookup_buckets = _partition_model_bucket_paths(
        inputs.partition_model_lookups, params.partition_lookup_keys
    )
    header = ["testing_unit_id", *[str(value) for value in range(1, params.bin_count + 1)]]
    components: dict[str, object] = {}
    model_outputs = dict(zip(expected_keys, outputs.partition_lookups, strict=True))
    for key, path in model_outputs.items():
        partition_id, model_id = key.split("\t", 1)
        model = next(model for model in enabled_models(params.hemizygote) if model.identifier == model_id)
        allowed = {
            unit for unit, partition in eligible_pairs
            if partition == partition_id
            and any(candidate.identifier == model_id for candidate in models_for_chromosome(
                next(row["CHROM"] for row in status_rows if row["testing_unit_id"] == unit and row["partition_id"] == partition),
                params.hemizygote,
            ))
        }
        _, inventory = write_gzip_tsv_atomic(path, header, lookup_rows(lookup_buckets[(partition_id, model_id)], allowed))
        components[f"lookup.{partition_id}.{model_id}"] = inventory

    if inputs.resolved_partition_intervals is None:
        raise ValueError("Partitioned bundle requires resolved ownership intervals")
    resolved = pl.read_parquet(inputs.resolved_partition_intervals)
    for partition_id, path in zip(params.partition_ids, outputs.resolved_partitions, strict=True):
        _, inventory = _write_partition_bed(
            path,
            resolved.filter(pl.col("partition_id") == partition_id)
            .select("CHROM", "START", "END", "testing_unit_id").iter_rows(),
        )
        components[f"resolved.{partition_id}"] = inventory
    for partition_id, path in zip(params.partition_ids, outputs.eligible_partitions, strict=True):
        _, inventory = _write_partition_bed(
            path,
            resolved.filter(
                (pl.col("partition_id") == partition_id)
                & pl.struct("testing_unit_id", "partition_id").is_in(
                    [{"testing_unit_id": unit, "partition_id": partition} for unit, partition in eligible_pairs]
                )
            ).select("CHROM", "START", "END", "testing_unit_id").iter_rows(),
        )
        components[f"eligible.{partition_id}"] = inventory
    status_header = list(status_rows[0]) if status_rows else ["testing_unit_id", "partition_id", "CHROM", "status", "failed_roles_json"]
    _, components["unit_partition_status"] = write_gzip_tsv_atomic(
        outputs.partition_statuses, status_header,
        ([row.get(name) for name in status_header] for row in status_rows),
    )
    overlap = pl.read_parquet(inputs.testing_unit_id_overlap).to_dicts() if inputs.testing_unit_id_overlap else []
    overlap_header = list(overlap[0]) if overlap else ["left_partition_id", "right_partition_id", "empty_set_status"]
    _, components["testing_unit_id_overlap"] = write_gzip_tsv_atomic(
        outputs.testing_unit_id_overlap, overlap_header,
        ([row.get(name) for name in overlap_header] for row in overlap),
    )
    calibration_count, components["calibration_regions"] = write_gzip_tsv_atomic(
        outputs.calibration_regions, ["#CHROM", "START", "END", "source"],
        parquet_rows([inputs.effective_mask], ["CHROM", "START", "END", "mask_state"]),
    )
    rejected_count, components["excluded_alleles"], components["candidate_audit"], zero_count, components["population_af_zero_imputations"], class_counts = write_finalization_allele_artifacts(
        inputs.allele_states, inputs.candidate_audit, outputs,
        include_variant_class_counts=params.variant_classes == ("snv", "indel"),
    )
    partition_metadata = []
    for index, partition_id in enumerate(params.partition_ids):
        rows = [row for row in status_rows if row["partition_id"] == partition_id]
        partition_metadata.append({
            "id": partition_id, "order": index,
            "source_identity": params.resource_identities[f"testing_unit_partition.{partition_id}"]["identity"],
            "resolved_component": f"resolved.{partition_id}",
            "eligible_component": f"eligible.{partition_id}",
            "lookups": {model.identifier: f"lookup.{partition_id}.{model.identifier}" for model in enabled_models(params.hemizygote)},
            "pair_counts": {"declared": len(rows), "eligible": sum(row["status"] == "eligible" for row in rows)},
        })
    metadata = {
        "schema_version": PARTITIONED_NULL_BUNDLE_SCHEMA_VERSION,
        "contract": PARTITIONED_NULL_BUNDLE_CONTRACT,
        "assembly": params.assembly, "contig_convention": CANONICAL_CONTIG_CONVENTION,
        "testing_unit_type": params.testing_unit_type, "af_info_tags": dict(params.af_tags),
        "binning": {"algorithm": FIXED_WIDTH_BINNING_ALGORITHM, "score_max": params.score_max, "bin_width": params.bin_width, "bin_count": params.bin_count, "maximum_in_final_bin": True},
        "model_assumptions": {"roles": list(model_roles(params.hemizygote)), "variant_classes": list(params.variant_classes), "unobserved_policy": params.policy, "mutation_conversion": params.conversion if params.policy == "mutation_rate" else None, "allele_events_independent": True, "par_policy": "excluded_before_event_assignment", "hemizygote_role": "hem_chrX" if params.hemizygote else None, "minimum_retained_fraction": params.minimum_retained_fraction, "retention_denominator": "resolved_supported_non_PAR_unit_bases_before_interval_exclusions", "x_chromosome": {"af_policy": params.x_af_policy, "mutation_rate_policy": params.x_mutation_rate_policy, "hemizygote_scope": "non_PAR_chrX_only"}},
        "partitions": partition_metadata,
        "overlap_policy": "earlier_partition_owns_overlapping_bases_v1",
        "indel_assignment_policy": "first_overlapping_resolved_partition_v1",
        "counts": {"eligible_pairs": len(eligible_pairs), "declared_pairs": len(status_rows), "excluded_alleles": rejected_count, "calibration_intervals": calibration_count, "population_af_zero_imputations": zero_count},
        "resource_identities": dict(params.resource_identities), "components": components,
        "validation_artifacts": dict(VALIDATION_ARTIFACT_FILENAMES),
    }
    if class_counts is not None:
        metadata["indel_universe"] = {"candidate_universe": params.indel_candidate_universe, "normalization": "vcf_minimal_left_aligned_v1", "reference_span": "zero_based_half_open_ref_span_v1", "membership": params.indel_membership, "mask_policy": params.indel_mask_policy, "max_reference_length": params.indel_max_reference_length, "max_alternate_length": params.indel_max_alternate_length}
        metadata["variant_class_counts"] = class_counts
    write_json_atomic(outputs.metadata, metadata)
