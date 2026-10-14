"""Window-local join of reference, score, and population observations."""

from __future__ import annotations

import polars as pl

from generation_bindings import (
    GenerationWindowInputs,
    GenerationWindowOutputs,
    GenerationWindowParameters,
)
from generation_windows import window_from_manifest
from generation_population import POPULATION_SCHEMA, population_rows
from generation_prep import intersect_intervals, merge_intervals
from generation_reference import ReferenceReader, membership_rows
from generation_scores import ALLELE_COLUMNS, SCORE_SCHEMA, score_rows
from generation_indels import (
    INDEL_REJECTION_SCHEMA,
    indel_candidate_memberships,
    indel_score_rows,
)
from io_contracts import write_parquet_atomic


def prepare_generation_window(
    inputs: GenerationWindowInputs,
    outputs: GenerationWindowOutputs,
    params: GenerationWindowParameters,
    window_id: str,
) -> None:
    """Materialize canonical states, memberships, and scores for one window."""
    state_schema = {**POPULATION_SCHEMA, "variant_class": pl.String}
    output_score_schema = {**SCORE_SCHEMA, "variant_class": pl.String}
    window = window_from_manifest(inputs.windows, window_id)
    halo_end = (
        window.end + params.indel_max_reference_length - 1
        if params.include_indels
        else window.end
    )
    eligible_units = (
        pl.scan_parquet(inputs.units)
        .filter(pl.col("status") == "eligible")
        .select("testing_unit_id")
        .collect(engine="streaming")
    )
    interval_frame = (
        pl.scan_parquet(inputs.unit_intervals)
        .filter(
            (pl.col("CHROM") == window.chromosome)
            & (pl.col("START") < halo_end)
            & (pl.col("END") > window.start)
        )
        .join(eligible_units.lazy(), on="testing_unit_id", how="inner")
        .select("testing_unit_id", "START", "END")
        .collect(engine="streaming")
    )
    mask_frame = (
        pl.scan_parquet(inputs.effective_mask)
        .filter(
            (pl.col("CHROM") == window.chromosome)
            & (pl.col("START") < halo_end)
            & (pl.col("END") > window.start)
        )
        .select("START", "END")
        .collect(engine="streaming")
    )
    intervals = interval_frame.filter(pl.col("START") < window.end).rows()
    mask = (
        mask_frame.filter(pl.col("START") < window.end)
        .select(
            pl.max_horizontal("START", pl.lit(window.start)).alias("START"),
            pl.min_horizontal("END", pl.lit(window.end)).alias("END"),
        )
        .rows()
    )
    reference = ReferenceReader(inputs.fasta, inputs.fasta_fai)
    try:
        snv_memberships = membership_rows(
            window, mask, intervals, reference
        ).with_columns(pl.lit("snv").alias("variant_class"))
        indel_scores = pl.DataFrame(schema=SCORE_SCHEMA)
        indel_memberships = pl.DataFrame(schema=snv_memberships.schema)
        indel_rejections = pl.DataFrame(schema=INDEL_REJECTION_SCHEMA)
        if params.include_indels:
            if inputs.score_indel is None:
                raise ValueError("Indel generation requires an indel score resource")
            all_intervals = interval_frame.rows()
            indel_mask = mask_frame.rows()
            catalog = indel_score_rows(
                inputs.score_indel,
                inputs.score_indel_index,
                window.chromosome,
                window.start,
                window.end,
                reference,
                max_reference_length=params.indel_max_reference_length,
                max_alternate_length=params.indel_max_alternate_length,
            )
            indel_memberships, indel_rejections = indel_candidate_memberships(
                catalog, indel_mask, all_intervals
            )
            indel_scores = catalog.select(*ALLELE_COLUMNS, "score")
    finally:
        reference.close()
    memberships = pl.concat(
        [snv_memberships, indel_memberships.select(*snv_memberships.columns)],
        how="vertical",
    )
    candidates = memberships.select(ALLELE_COLUMNS).unique(maintain_order=True)
    candidate_classes = memberships.select(*ALLELE_COLUMNS, "variant_class").unique(
        maintain_order=True
    )
    if candidates.is_empty():
        # Keep the same artifact contract without opening either large resource.
        write_parquet_atomic(
            pl.DataFrame(schema=state_schema), outputs.allele_states, params.compression
        )
        write_parquet_atomic(memberships, outputs.allele_membership, params.compression)
        write_parquet_atomic(
            pl.DataFrame(schema=output_score_schema),
            outputs.allele_scores,
            params.compression,
        )
        write_parquet_atomic(
            indel_rejections, outputs.candidate_audit, params.compression
        )
        return
    candidate_intervals = intersect_intervals(
        [(start, end) for _, start, end in intervals], mask
    )
    # Indel membership can be earned by a span that overlaps a unit even when
    # its VCF anchor is just outside that unit's callable bases.  Retain each
    # anchor as a one-base query region so the population join sees it too.
    indel_anchors = [
        (position - 1, position)
        for position in candidate_classes.filter(pl.col("variant_class") == "indel")[
            "POS"
        ]
        .unique()
        .to_list()
    ]
    candidate_intervals = merge_intervals(candidate_intervals + indel_anchors)
    scores = score_rows(
        inputs.score_snv,
        inputs.score_snv_index,
        window.chromosome,
        window.start,
        window.end,
        intervals=candidate_intervals,
    )
    scores = pl.concat([scores, indel_scores], how="vertical")
    matched_score_keys = candidates.join(
        scores.select(ALLELE_COLUMNS), on=ALLELE_COLUMNS, how="inner"
    )
    if matched_score_keys.select(ALLELE_COLUMNS).is_duplicated().any():
        raise ValueError(
            f"Score resource has duplicate calibration alleles in {window.identifier}"
        )
    scored = candidates.join(scores, on=ALLELE_COLUMNS, how="inner").join(
        candidate_classes, on=ALLELE_COLUMNS, how="inner"
    )
    if scored.height != candidates.height:
        missing = candidates.join(
            scores.select(ALLELE_COLUMNS), on=ALLELE_COLUMNS, how="anti"
        )
        examples = ", ".join(
            f"{row[0]}:{row[1]}:{row[2]}:{row[3]}"
            for row in missing.head(10).iter_rows()
        )
        raise ValueError(
            f"Score resource is missing {missing.height} calibration alleles in {window.identifier}: {examples}"
        )
    panel = population_rows(
        inputs.population_vcf,
        inputs.population_vcf_index,
        window.chromosome,
        window.start,
        window.end,
        dict(params.af_tags),
        params.x_af_policy,
        params.hemizygote,
        intervals=candidate_intervals,
    )
    if panel.select(ALLELE_COLUMNS).is_duplicated().any():
        duplicates = panel.filter(pl.struct(ALLELE_COLUMNS).is_duplicated())
        examples = ", ".join(
            f"{chromosome}:{position}:{reference}:{alternate}"
            for chromosome, position, reference, alternate in duplicates.head(10)
            .select(ALLELE_COLUMNS)
            .iter_rows()
        )
        raise ValueError(
            f"Population VCF has duplicate allele states in {window.identifier}: {examples}"
        )
    states = (
        candidates.join(panel, on=ALLELE_COLUMNS, how="left")
        .with_columns(
            pl.when(pl.col("allele_state").is_null())
            .then(pl.lit("absent"))
            .otherwise(pl.col("allele_state"))
            .alias("allele_state"),
            pl.when(pl.col("state_reason").is_null())
            .then(pl.lit("not_present_in_population_panel"))
            .otherwise(pl.col("state_reason"))
            .alias("state_reason"),
        )
        .join(candidate_classes, on=ALLELE_COLUMNS, how="inner")
        .select(
            *ALLELE_COLUMNS,
            "allele_state",
            "af",
            "af_overall",
            "af_female",
            "af_male",
            "state_reason",
            "af_zero_imputed_tags",
            "variant_class",
        )
    )
    invalid = states.filter(pl.col("allele_state") == "invalid")
    if invalid.height:
        examples = "; ".join(
            f"{chrom}:{position}:{reference}:{alternate}: {reason}"
            for chrom, position, reference, alternate, reason in invalid.head(10)
            .select(*ALLELE_COLUMNS, "state_reason")
            .iter_rows()
        )
        raise ValueError(
            f"Population VCF has invalid selected AF for {invalid.height} calibration allele(s) in {window.identifier}: {examples}"
        )
    write_parquet_atomic(states, outputs.allele_states, params.compression)
    write_parquet_atomic(memberships, outputs.allele_membership, params.compression)
    write_parquet_atomic(scored, outputs.allele_scores, params.compression)
    write_parquet_atomic(indel_rejections, outputs.candidate_audit, params.compression)
