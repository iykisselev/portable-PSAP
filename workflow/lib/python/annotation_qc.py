#!/usr/bin/env python3
"""Summarize whether cohort variants received CADD and testing-unit annotations."""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

from io_contracts import is_missing_token


REQUIRED_COLUMNS = ("CHROM", "POS", "REF", "ALT", "CADD_PHRED", "unit")


def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref.upper() in "ACGT" and alt.upper() in "ACGT"


def has_value(value: str) -> bool:
    return not is_missing_token(value)


def summarize_annotations(path: str | Path, include_indels: bool) -> dict:
    errors: list[str] = []
    counts = {
        "variant_rows": 0,
        "snv_rows": 0,
        "indel_rows": 0,
        "cadd_annotated_rows": 0,
        "snv_cadd_annotated_rows": 0,
        "indel_cadd_annotated_rows": 0,
        "testing_unit_annotated_rows": 0,
        "fully_annotated_rows": 0,
        "invalid_cadd_rows": 0,
    }

    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing_columns = [
                column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])
            ]
            if missing_columns:
                errors.append(
                    "Annotation table is missing required column(s): "
                    + ", ".join(missing_columns)
                )
            else:
                for row in reader:
                    counts["variant_rows"] += 1
                    snv = is_snv(row["REF"].strip(), row["ALT"].strip())
                    variant_key = "snv" if snv else "indel"
                    counts[f"{variant_key}_rows"] += 1

                    cadd = has_value(row["CADD_PHRED"])
                    if cadd:
                        try:
                            cadd = math.isfinite(float(row["CADD_PHRED"]))
                        except ValueError:
                            cadd = False
                        if not cadd:
                            counts["invalid_cadd_rows"] += 1
                    unit = has_value(row["unit"])

                    if cadd:
                        counts["cadd_annotated_rows"] += 1
                        counts[f"{variant_key}_cadd_annotated_rows"] += 1
                    if unit:
                        counts["testing_unit_annotated_rows"] += 1
                    if cadd and unit:
                        counts["fully_annotated_rows"] += 1
    except (OSError, EOFError) as exc:
        errors.append(f"Could not read annotation table {path}: {exc}")

    if counts["invalid_cadd_rows"]:
        errors.append(
            f"Annotation table contains {counts['invalid_cadd_rows']} non-numeric "
            "or non-finite CADD_PHRED value(s)"
        )
    if counts["snv_rows"] and not counts["snv_cadd_annotated_rows"]:
        errors.append(
            "CADD SNV annotation matched zero SNV records; check the CADD assembly "
            "and chromosome naming convention"
        )
    if (
        include_indels
        and counts["indel_rows"]
        and not counts["indel_cadd_annotated_rows"]
    ):
        errors.append(
            "CADD indel annotation matched zero indel records; check the CADD "
            "assembly and chromosome naming convention"
        )
    if counts["variant_rows"] and not counts["testing_unit_annotated_rows"]:
        errors.append(
            "Testing-unit annotation matched zero variant records; check the BED "
            "assembly, coordinates, and chromosome naming convention"
        )

    return {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": [],
        "summary": counts,
    }
