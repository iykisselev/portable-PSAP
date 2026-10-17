"""Audit testing-unit identifier overlap among configured partitions."""

from __future__ import annotations

from typing import Mapping


def partition_id_overlap_rows(
    partitions: list[tuple[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    """Return one ordered pair record with explicit empty-set semantics."""
    rows: list[dict[str, object]] = []
    for left_index, (left_id, left_units) in enumerate(partitions):
        left = set(left_units)
        for right_id, right_units in partitions[left_index + 1 :]:
            right = set(right_units)
            intersection = len(left & right)
            union = len(left | right)
            minimum = min(len(left), len(right))
            rows.append(
                {
                    "left_partition_id": left_id,
                    "right_partition_id": right_id,
                    "left_unique_id_count": len(left),
                    "right_unique_id_count": len(right),
                    "intersection_count": intersection,
                    "union_count": union,
                    "left_overlap_fraction": intersection / len(left) if left else None,
                    "right_overlap_fraction": intersection / len(right) if right else None,
                    "minimum_overlap_fraction": intersection / minimum if minimum else None,
                    "empty_set_status": (
                        "both_empty"
                        if not left and not right
                        else "left_empty"
                        if not left
                        else "right_empty"
                        if not right
                        else "nonempty"
                    ),
                    "low_overlap_warning": bool(minimum and intersection / minimum < 0.5),
                }
            )
    return rows
