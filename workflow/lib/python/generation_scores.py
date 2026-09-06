"""Partition-local score and mutation-rate readers."""

from __future__ import annotations

import math
from pathlib import Path

import polars as pl

try:
    from generation_population import indexed_lines
    from genomic_contracts import allele_key
except ModuleNotFoundError:
    from .generation_population import indexed_lines
    from .genomic_contracts import allele_key


ALLELE_COLUMNS = ["CHROM", "POS", "REF", "ALT"]


def score_rows(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
) -> pl.DataFrame:
    """Read CADD-layout SNV scores from one genomic window."""
    records: list[dict[str, object]] = []
    for raw in indexed_lines(path, index, chromosome, start, end):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise ValueError("Score resource needs at least six CADD columns")
        key = allele_key(fields[0], int(fields[1]), fields[2], fields[3])
        if key[0] != chromosome or not start < key[1] <= end or len(key[2]) != 1 or len(key[3]) != 1:
            continue
        score = float(fields[5])
        if not math.isfinite(score):
            raise ValueError(f"Score resource has a non-finite score for {key!r}")
        records.append({"CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3], "score": score})
    return pl.DataFrame(
        records,
        schema={"CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String, "score": pl.Float64},
        strict=False,
    )


def mutation_rate_rows(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
) -> pl.DataFrame:
    """Read per-allele mutation rates from one requested genomic window."""
    records: list[dict[str, object]] = []
    for raw in indexed_lines(path, index, chromosome, start, end):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 5:
            raise ValueError("Mutation-rate resource needs CHROM, POS, REF, ALT, and mu")
        key = allele_key(fields[0], int(fields[1]), fields[2], fields[3])
        if key[0] != chromosome or not start < key[1] <= end:
            continue
        rate = float(fields[4])
        if not math.isfinite(rate) or not 0 <= rate <= 1:
            raise ValueError(f"Mutation-rate resource has invalid mu for {key!r}")
        records.append({"CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3], "mu": rate})
    return pl.DataFrame(
        records,
        schema={"CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String, "mu": pl.Float64},
        strict=False,
    )
