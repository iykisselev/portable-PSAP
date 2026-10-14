"""Window-manifest adapters for bounded null-distribution generation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from generation_manifests import read_window_manifest
from genomic_contracts import read_fai


@dataclass(frozen=True)
class GenomeWindow:
    """One half-open genomic region with a stable output identifier."""

    identifier: str
    chromosome: str
    start: int
    end: int
    ordinal: int


def read_fasta_index(path: str | Path) -> dict[str, tuple[str, int]]:
    """Compatibility view of the shared FAI records for windowing."""
    return {
        chromosome: (record.name, record.length)
        for chromosome, record in read_fai(path).items()
    }


def build_windows(
    fasta_index: str | Path,
    testing_unit_intervals: list[tuple[str, int, int]],
    window_size: int,
    *,
    candidate_owner_intervals: list[tuple[str, int, int]] | None = None,
) -> list[GenomeWindow]:
    """Return fixed-size windows that intersect eligible candidate intervals.

    ``candidate_owner_intervals`` is an optional indel anchor halo.  It is
    separate from the ordinary callable-unit intervals so downstream code can
    preserve exact span membership while still scheduling a preceding owner
    window for a boundary-crossing deletion.
    """
    if window_size <= 0:
        raise ValueError("generation.execution.partition_size_bp must be positive")
    fasta_contigs = read_fasta_index(fasta_index)
    candidate_intervals = list(testing_unit_intervals)
    if candidate_owner_intervals is not None:
        candidate_intervals.extend(candidate_owner_intervals)
    chromosomes = {chromosome for chromosome, _, _ in candidate_intervals}
    missing = sorted(chromosomes - set(fasta_contigs))
    if missing:
        raise ValueError(
            "Testing-unit contigs are absent from the FASTA index: "
            + ", ".join(missing)
        )

    selected_windows: dict[str, set[int]] = {}
    for chromosome, interval_start, interval_end in candidate_intervals:
        _, contig_length = fasta_contigs[chromosome]
        if interval_start < 0 or interval_end <= interval_start:
            raise ValueError(f"Testing-unit interval on {chromosome} is invalid")
        if interval_end > contig_length:
            raise ValueError(
                f"Testing-unit interval on {chromosome} extends beyond the FASTA contig"
            )
        first_window = interval_start // window_size
        last_window = (interval_end - 1) // window_size
        selected_windows.setdefault(chromosome, set()).update(
            range(first_window, last_window + 1)
        )

    windows: list[GenomeWindow] = []
    ordinal = 0
    for chromosome in sorted(selected_windows):
        _, length = fasta_contigs[chromosome]
        for window in sorted(selected_windows[chromosome]):
            start = window * window_size
            end = min(start + window_size, length)
            windows.append(
                GenomeWindow(
                    identifier=f"{chromosome}_{start:012d}_{end:012d}",
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    ordinal=ordinal,
                )
            )
            ordinal += 1
    return windows


def indel_candidate_owner_intervals(
    testing_unit_intervals: list[tuple[str, str, int, int]],
    max_reference_length: int,
) -> list[tuple[str, str, int, int]]:
    """Expand unit intervals leftward for indel anchor-owner discovery.

    An indel with a REF of at most ``max_reference_length`` can overlap a unit
    while its anchor is up to ``max_reference_length - 1`` bases upstream.
    The returned intervals are not membership intervals: exact span overlap is
    still checked by :func:`generation_indels.indel_candidate_memberships`.
    """
    if max_reference_length < 1:
        raise ValueError("Indel maximum reference length must be positive")
    return [
        (unit_id, chromosome, max(0, start - max_reference_length + 1), end)
        for unit_id, chromosome, start, end in testing_unit_intervals
    ]


def read_genome_windows(path: str | Path) -> list[GenomeWindow]:
    """Load validated scope windows into the window-local record type."""
    return [
        GenomeWindow(
            identifier=entry.identifier,
            chromosome=entry.chromosome,
            start=entry.start,
            end=entry.end,
            ordinal=entry.ordinal,
        )
        for entry in read_window_manifest(path)
    ]


def window_from_manifest(path: str | Path, identifier: str) -> GenomeWindow:
    """Resolve one requested wildcard without reconstructing coordinates from text."""
    for window in read_genome_windows(path):
        if window.identifier == identifier:
            return window
    raise ValueError(f"Unknown generation window {identifier!r}")
