"""Streaming assembly of generated null-distribution bundle components."""

from __future__ import annotations

import csv
from collections import defaultdict
from itertools import chain
import math
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
    INDEL_GENERATED_BUNDLE_CONTRACT,
    INDEL_GENERATED_BUNDLE_SCHEMA_VERSION,
    VALIDATION_ARTIFACT_FILENAMES,
)
from io_contracts import write_gzip_tsv_atomic, write_json_atomic


def applicable_roles(unit_row, hemizygote: bool) -> tuple[str, ...]:
    """Return roles needed to score a testing unit under the X-only HEM contract."""
    return ("het", "hom", "hem_chrX") if hemizygote and unit_row["CHROM"] == "X" else ("het", "hom")


def classify_units(units, diagnostics, hemizygote: bool = False):
    """Keep units with positive mass in every role applicable to their chromosome."""
    unit_status = {row["testing_unit_id"]: row["status"] for row in units}
    if len(unit_status) != len(units):
        raise ValueError("Provisional unit table contains duplicate testing-unit IDs")
    by_unit = defaultdict(dict)
    for row in diagnostics:
        unit, role, status = row["testing_unit_id"], row["role"], row["status"]
        if unit not in unit_status:
            raise ValueError(f"Distribution diagnostics contain unknown unit {unit!r}")
        if role not in {"het", "hom", "hem_chrX"}:
            raise ValueError(f"Distribution diagnostics contain unsupported role {role!r}")
        if role in by_unit[unit]:
            raise ValueError(f"Duplicate distribution diagnostic for {unit}/{role}")
        try:
            event_mass = float(row["event_mass"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid event mass for {unit}/{role}") from exc
        if status == "ok" and (not math.isfinite(event_mass) or not 0 < event_mass <= 1):
            raise ValueError(f"Invalid positive event mass for {unit}/{role}")
        if status == "zero_event_mass" and event_mass != 0:
            raise ValueError(f"Zero-event diagnostic has nonzero mass for {unit}/{role}")
        if status == "numerical_failure":
            raise ValueError(f"Numerical distribution failure for {unit}/{role}")
        if status not in {"ok", "zero_event_mass", "numerical_failure"}:
            raise ValueError(f"Unsupported distribution status {status!r} for {unit}/{role}")
        by_unit[unit][role] = status
    eligible, excluded = [], []
    for row in units:
        unit = row["testing_unit_id"]
        if row["status"] != "eligible":
            excluded.append((unit, row["status"]))
            continue
        roles = applicable_roles(row, hemizygote)
        missing = [role for role in roles if role not in by_unit[unit]]
        if missing:
            raise ValueError(f"Missing distribution diagnostic for {unit}: {', '.join(missing)}")
        zero = [role for role in roles if by_unit[unit][role] == "zero_event_mass"]
        if zero:
            excluded.append((unit, "zero_event_mass_" + "_".join(zero)))
        else:
            eligible.append(unit)
    return eligible, excluded


def diagnostics_rows(paths):
    """Read compact bucket diagnostics; their size is bounded by testing units."""
    for path in paths:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


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
        yield from pl.scan_parquet(paths).select(columns).collect(engine="streaming").iter_rows()
        return
    for path in paths:
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(columns=columns):
            values = batch.to_pydict()
            yield from zip(*(values[column] for column in columns))


def masked_bed_rows(path: str | Path, allowed_units):
    """Stream canonical masked intervals, never the raw pre-PAR testing BED."""
    for chromosome, start, end, unit in parquet_rows([path], ["CHROM", "START", "END", "testing_unit_id"]):
        if unit in allowed_units:
            yield chromosome, start, end, unit


def rejected_rows(paths):
    """Expose rejected states without loading all partition state tables."""
    for chromosome, position, reference, alternate, state, reason in parquet_rows(
        paths, ["CHROM", "POS", "REF", "ALT", "allele_state", "state_reason"]
    ):
        if state == "rejected":
            yield chromosome, position, reference, alternate, reason


def candidate_audit_rows(paths):
    """Expose structural indel candidate rejections for the public audit."""
    for chromosome, position, reference, alternate, reason in parquet_rows(
        paths, ["CHROM", "POS", "REF", "ALT", "rejection_reason"]
    ):
        yield chromosome, position, reference, alternate, reason


def _empty_variant_class_counts() -> dict[str, dict[str, int]]:
    return {
        variant_class: {name: 0 for name in ("retained", "rejected", "observed", "absent")}
        for variant_class in ("snv", "indel")
    }


def _update_variant_class_counts(
    counts: dict[str, dict[str, int]], variant_class: str, state: str, amount: int = 1,
) -> None:
    """Add one grouped population-state count to the v3 metadata counts."""
    if variant_class not in counts:
        raise ValueError(f"Unsupported variant class in allele states: {variant_class!r}")
    if state in {"observed", "absent"}:
        counts[variant_class]["retained"] += amount
        counts[variant_class][state] += amount
    elif state == "rejected":
        counts[variant_class]["rejected"] += amount
    else:
        raise ValueError(f"Unsupported allele state in bundle finalization: {state!r}")


def _variant_class_counts_from_frame(count_frame: pl.DataFrame) -> dict[str, dict[str, int]]:
    """Build public v3 counts from a compact grouped state-count frame."""
    counts = _empty_variant_class_counts()
    for variant_class, state, amount in count_frame.iter_rows():
        _update_variant_class_counts(counts, variant_class, state, amount)
    return counts


def variant_class_counts(state_paths, audit_paths) -> dict[str, dict[str, int]]:
    """Count retained, population-rejected, and structural candidates by class."""
    count_frame = (
        pl.scan_parquet(state_paths)
        .group_by("variant_class", "allele_state")
        .len()
        .collect(engine="streaming")
    )
    counts = _variant_class_counts_from_frame(count_frame)
    for _ in candidate_audit_rows(audit_paths):
        counts["indel"]["rejected"] += 1
    return counts


def zero_imputed_af_rows(paths):
    """Expose calibration alleles whose omitted population AF tag became zero."""
    for chromosome, position, reference, alternate, tags in parquet_rows(
        paths, ["CHROM", "POS", "REF", "ALT", "af_zero_imputed_tags"]
    ):
        if tags:
            yield chromosome, position, reference, alternate, tags


def collect_finalization_allele_rows(
    state_paths,
    audit_paths,
    *,
    include_variant_class_counts: bool,
) -> tuple[
    list[tuple[object, ...]],
    list[tuple[object, ...]],
    list[tuple[object, ...]],
    dict[str, dict[str, int]] | None,
]:
    """Collect sparse public rows and grouped counts with native Parquet scans.

    Candidate-state shards can contain hundreds of millions of rows, while
    public rejection and zero-imputation artifacts are normally tiny. Keep the
    candidate scan in Polars and materialize only those exceptional rows.
    State and audit lists remain separate so excluded-allele ordering stays
    state-first followed by structural indel rejections.
    """
    states = pl.scan_parquet(state_paths)
    scans = [
        states.filter(pl.col("allele_state") == "rejected").select(
            "CHROM", "POS", "REF", "ALT", "state_reason"
        ),
        states.filter(
            pl.col("af_zero_imputed_tags").is_not_null()
            & (pl.col("af_zero_imputed_tags") != "")
        ).select("CHROM", "POS", "REF", "ALT", "af_zero_imputed_tags"),
    ]
    if include_variant_class_counts:
        scans.append(states.group_by("variant_class", "allele_state").len())
    collected = pl.collect_all(scans, engine="streaming")
    rejected = list(collected[0].iter_rows())
    zero_imputed = list(collected[1].iter_rows())
    audit = list(candidate_audit_rows(audit_paths))
    counts = (
        _variant_class_counts_from_frame(collected[2])
        if include_variant_class_counts
        else None
    )
    if counts is not None:
        counts["indel"]["rejected"] += len(audit)
    return rejected, zero_imputed, audit, counts


def finalize_generated_bundle(
    inputs: BundleFinalizationInputs,
    outputs: BundleFinalizationOutputs,
    params: BundleFinalizationParameters,
) -> None:
    """Build a certified public bundle from bounded generation artifacts."""
    combined_indel_bundle = params.variant_classes == ("snv", "indel")
    if params.variant_classes not in {("snv",), ("snv", "indel")}:
        raise ValueError(f"Unsupported generated variant classes: {params.variant_classes!r}")
    units = pl.read_parquet(inputs.units).to_dicts()
    eligible, excluded = classify_units(units, list(diagnostics_rows(inputs.diagnostics)), hemizygote=params.hemizygote)
    eligible_set = set(eligible)
    eligible_x_units = {row["testing_unit_id"] for row in units if row["CHROM"] == "X" and row["testing_unit_id"] in eligible_set}
    if params.hemizygote and not eligible_x_units:
        raise ValueError("HEM generation requested, but no chromosome-X testing unit survived the HET/HOM/HEM eligibility contract")
    header = ["testing_unit_id", *[str(value) for value in range(1, params.bin_count + 1)]]
    _, lookup_het = write_gzip_tsv_atomic(outputs.lookup_het, header, lookup_rows(inputs.lookup_het, eligible_set))
    _, lookup_hom = write_gzip_tsv_atomic(outputs.lookup_hom, header, lookup_rows(inputs.lookup_hom, eligible_set))
    lookup_hem = None
    if params.hemizygote:
        hem_count, lookup_hem = write_gzip_tsv_atomic(outputs.lookup_hem, header, lookup_rows(inputs.lookup_hem, eligible_x_units))
        if hem_count != len(eligible_x_units):
            raise ValueError("HEM lookup rows do not exactly cover surviving chromosome-X testing units")
    _, excluded_units = write_gzip_tsv_atomic(outputs.excluded_units, ["testing_unit_id", "exclusion_reason"], excluded)
    _, eligible_bed = write_gzip_tsv_atomic(outputs.eligible_bed, ["#CHROM", "START", "END", "testing_unit_id"], masked_bed_rows(inputs.masked_unit_intervals, eligible_set))
    calibration_count, calibration_regions = write_gzip_tsv_atomic(
        outputs.calibration_regions, ["#CHROM", "START", "END", "source"], parquet_rows([inputs.effective_mask], ["CHROM", "START", "END", "mask_state"])
    )
    rejected_state_rows, zero_imputed_rows, audit_rows, class_counts = collect_finalization_allele_rows(
        inputs.allele_states,
        inputs.candidate_audit,
        include_variant_class_counts=combined_indel_bundle,
    )
    rejected_count, excluded_alleles = write_gzip_tsv_atomic(
        outputs.excluded_alleles, ["#CHROM", "POS", "REF", "ALT", "exclusion_reason"],
        chain(rejected_state_rows, audit_rows)
    )
    _, candidate_audit = write_gzip_tsv_atomic(
        outputs.candidate_audit, ["#CHROM", "POS", "REF", "ALT", "rejection_reason"],
        audit_rows,
    )
    zero_imputation_count, population_af_zero_imputations = write_gzip_tsv_atomic(
        outputs.population_af_zero_imputations,
        ["#CHROM", "POS", "REF", "ALT", "imputed_info_tags"],
        zero_imputed_rows,
    )
    components = {
        "lookup_het": lookup_het, "lookup_hom": lookup_hom,
        "eligible_testing_units": eligible_bed, "calibration_regions": calibration_regions,
        "excluded_alleles": excluded_alleles,
        "candidate_audit": candidate_audit,
        "population_af_zero_imputations": population_af_zero_imputations,
        "excluded_units": excluded_units,
    }
    if params.hemizygote:
        assert lookup_hem is not None
        components["lookup_hem_chrX"] = lookup_hem
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
        "binning": {"algorithm": FIXED_WIDTH_BINNING_ALGORITHM, "score_max": params.score_max, "bin_width": params.bin_width, "bin_count": params.bin_count, "maximum_in_final_bin": True},
        "model_assumptions": {
            "minimum_retained_fraction": params.minimum_retained_fraction,
            "retention_denominator": "supported_non_PAR_unit_bases_before_interval_exclusions",
            "roles": ["het", "hom"] + (["hem_chrX"] if params.hemizygote else []),
            "variant_classes": list(params.variant_classes), "unobserved_policy": params.policy,
            "mutation_conversion": params.conversion if params.policy == "mutation_rate" else None,
            "allele_events_independent": True, "par_policy": "excluded_before_event_assignment",
            "hemizygote_role": "hem_chrX" if params.hemizygote else None,
            "x_chromosome": {"af_policy": params.x_af_policy, "mutation_rate_policy": params.x_mutation_rate_policy, "hemizygote_scope": "non_PAR_chrX_only"},
        },
        "counts": {
            "eligible_units": len(eligible), "excluded_units": len(excluded),
            "excluded_alleles": rejected_count, "calibration_intervals": calibration_count,
            "population_af_zero_imputations": zero_imputation_count,
        },
        "resource_identities": dict(params.resource_identities), "components": components,
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
