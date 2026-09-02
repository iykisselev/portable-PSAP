"""Compatibility exports for the Stage 6 execution-policy module."""

from __future__ import annotations

from typing import Iterable

from execution_policy import input_size_mb, resource_mem_mb, resource_runtime_minutes


def scaled_mem_mb(
    input_files: Iterable[object], attempt: int, baseline_mb: int, mb_per_gb: int = 0
) -> int:
    """Preserve the legacy scalar interface for external rule adapters."""
    from execution_policy import scaled_value

    return scaled_value(baseline_mb, mb_per_gb, input_files, attempt)


def scaled_runtime_minutes(
    input_files: Iterable[object],
    attempt: int,
    baseline_minutes: int,
    minutes_per_gb: int = 0,
) -> int:
    """Preserve the legacy scalar interface for external rule adapters."""
    from execution_policy import scaled_value

    return scaled_value(baseline_minutes, minutes_per_gb, input_files, attempt)
