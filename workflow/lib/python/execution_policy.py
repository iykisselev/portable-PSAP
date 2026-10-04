"""Scheduler-independent, named execution policies for workflow rules.

Rules retain explicit scheduler and software declarations. They select only a
logical policy name here, so baseline resources and size-based retry scaling
are tested in one place rather than duplicated as numeric multipliers.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from types import MappingProxyType
from typing import Iterable, Mapping


@dataclass(frozen=True)
class ExecutionPolicy:
    """Baseline memory/runtime plus whole-GiB input scaling coefficients."""

    baseline_mem_mb: int
    baseline_runtime_minutes: int
    mem_mb_per_gib: int = 0
    runtime_minutes_per_gib: int = 0


EXECUTION_POLICIES: Mapping[str, ExecutionPolicy] = MappingProxyType(
    {
        "stage_default": ExecutionPolicy(4096, 180),
        "stage_vcf": ExecutionPolicy(4096, 180),
        "stage_calibration_regions": ExecutionPolicy(2048, 60),
        "stage_testing_units": ExecutionPolicy(2048, 60),
        "stage_cadd_snv": ExecutionPolicy(2048, 360, 256, 20),
        "stage_cadd_indel": ExecutionPolicy(2048, 360, 256, 20),
        "annotate_variants": ExecutionPolicy(32000, 720, 512, 30),
        "cohort_aggregate": ExecutionPolicy(16000, 360, 256, 20),
    }
)


def input_size_mb(input_files: Iterable[object]) -> float:
    """Return declared local-input bytes in MiB, ignoring unavailable paths."""
    size_mb = 0.0
    for item in input_files:
        try:
            size_mb += Path(str(item)).stat().st_size / (1024 * 1024)
        except OSError:
            pass
    return size_mb


def execution_policy(name: str) -> ExecutionPolicy:
    """Resolve a named policy, using the shared staging fallback when needed."""
    try:
        return EXECUTION_POLICIES[name]
    except KeyError:
        if name.startswith("stage_"):
            return EXECUTION_POLICIES["stage_default"]
        raise ValueError(f"Unknown execution policy: {name!r}") from None


def scaled_value(baseline: int, per_gib: int, input_files: Iterable[object], attempt: int) -> int:
    """Apply whole-GiB input scaling and Snakemake retry scaling."""
    return int(
        (baseline + math.ceil(input_size_mb(input_files) / 1024) * per_gib)
        * attempt
    )


def resource_mem_mb(name: str, input_files: Iterable[object], attempt: int) -> int:
    """Return memory requested for one logical execution policy."""
    policy = execution_policy(name)
    return scaled_value(
        policy.baseline_mem_mb,
        policy.mem_mb_per_gib,
        input_files,
        attempt,
    )


def resource_runtime_minutes(name: str, input_files: Iterable[object], attempt: int) -> int:
    """Return wall time requested for one logical execution policy."""
    policy = execution_policy(name)
    return scaled_value(
        policy.baseline_runtime_minutes,
        policy.runtime_minutes_per_gib,
        input_files,
        attempt,
    )
