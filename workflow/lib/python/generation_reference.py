"""Reference sequence and testing-unit membership construction."""

from __future__ import annotations

from bisect import bisect_right
from pathlib import Path

import numpy as np
import polars as pl

from generation_windows import GenomeWindow
from generation_prep import merge_intervals
from genomic_contracts import read_fai

try:
    import pysam
except ImportError:
    pysam = None


MEMBERSHIP_SCHEMA = {
    "testing_unit_id": pl.String,
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
}
_BASE_CODES = np.full(256, -1, dtype=np.int8)
_BASE_CODES[[ord("A"), ord("C"), ord("G"), ord("T")]] = np.arange(4, dtype=np.int8)
_BASES = np.array(["A", "C", "G", "T"])
# Rows are ordered exactly as ``for alternate in 'ACGT'`` after excluding REF.
_ALTERNATES = np.array(
    [
        ["C", "G", "T"],
        ["A", "G", "T"],
        ["A", "C", "T"],
        ["A", "C", "G"],
    ]
)


class ReferenceReader:
    """Random-access FASTA reader using pysam or uncompressed FAI offsets."""

    def __init__(self, fasta: str | Path, fai: str | Path) -> None:
        self.fasta = Path(fasta)
        self.fai = Path(fai)
        self.contigs = read_fai(self.fai)
        self.reader = (
            pysam.FastaFile(str(self.fasta), filepath_index=str(self.fai))
            if pysam is not None
            else None
        )
        if self.reader is None:
            with self.fasta.open("rb") as handle:
                if handle.read(2) == b"\x1f\x8b":
                    raise RuntimeError("Compressed FASTA requires pysam")

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        """Fetch one bounded half-open interval with canonical chromosome names."""
        record = self.contigs[chromosome]
        if start < 0 or end > record.length or start > end:
            raise ValueError(
                f"Reference request is outside {chromosome!r}: {start}-{end}"
            )
        if self.reader is not None:
            sequence = self.reader.fetch(record.name, start, end).upper()
            if len(sequence) != end - start:
                raise ValueError(
                    f"FASTA is truncated or its index is stale: {self.fasta}"
                )
            return sequence

        sequence: list[bytes] = []
        position = start
        with self.fasta.open("rb") as handle:
            while position < end:
                within_line = position % record.line_bases
                take = min(end - position, record.line_bases - within_line)
                handle.seek(
                    record.offset
                    + (position // record.line_bases) * record.line_bytes
                    + within_line
                )
                chunk = handle.read(take)
                if len(chunk) != take:
                    raise ValueError(
                        f"FASTA is truncated or its index is stale: {self.fasta}"
                    )
                sequence.append(chunk)
                position += take
        return b"".join(sequence).decode("ascii").upper()

    def close(self) -> None:
        """Close the indexed reader when pysam owns an open file handle."""
        if self.reader is not None:
            self.reader.close()


def membership_rows(
    window: GenomeWindow,
    mask: list[tuple[int, int]],
    intervals: list[tuple[str, int, int]],
    reference: ReferenceReader,
) -> pl.DataFrame:
    """Generate bounded candidate SNVs and distinct allele-to-unit memberships."""
    frames: list[pl.DataFrame] = []
    if any(mask[index - 1][1] > mask[index][0] for index in range(1, len(mask))):
        mask = merge_intervals(mask)
    mask_starts = [start for start, _ in mask]
    # Canonical scope artifacts have disjoint masks and merged intervals for
    # each unit, so no two chunks can emit the same membership.  Retaining a
    # full-frame unique operation in that normal case would create an avoidable
    # hash table at the largest point in a window job.  Keep the public
    # distinct-row behavior for direct callers with overlapping inputs.
    greatest_end_by_unit: dict[str, int] = {}
    needs_deduplication = False
    for unit, start, end in intervals:
        if start < greatest_end_by_unit.get(unit, start):
            needs_deduplication = True
            break
        greatest_end_by_unit[unit] = max(end, greatest_end_by_unit.get(unit, end))
    for unit, start, end in intervals:
        mask_index = max(0, bisect_right(mask_starts, start) - 1)
        while mask_index < len(mask) and mask[mask_index][0] < end:
            mask_start, mask_end = mask[mask_index]
            overlap_start, overlap_end = max(start, mask_start), min(end, mask_end)
            mask_index += 1
            if overlap_start >= overlap_end:
                continue
            sequence = reference.fetch(window.chromosome, overlap_start, overlap_end)
            encoded = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
            base_codes = _BASE_CODES[encoded]
            valid = base_codes >= 0
            if not valid.any():
                continue
            valid_codes = base_codes[valid]
            positions = np.arange(overlap_start, overlap_end, dtype=np.int64)[valid]
            frames.append(
                pl.DataFrame(
                    {
                        "POS": np.repeat(positions + 1, 3),
                        "REF": np.repeat(_BASES[valid_codes], 3),
                        "ALT": _ALTERNATES[valid_codes].reshape(-1),
                    },
                    schema={"POS": pl.Int64, "REF": pl.String, "ALT": pl.String},
                    strict=False,
                )
                .with_columns(
                    pl.lit(unit).alias("testing_unit_id"),
                    pl.lit(window.chromosome).alias("CHROM"),
                )
                .select(*MEMBERSHIP_SCHEMA)
            )
    if not frames:
        return pl.DataFrame(schema=MEMBERSHIP_SCHEMA)
    memberships = pl.concat(frames, how="vertical")
    return (
        memberships.unique(maintain_order=True) if needs_deduplication else memberships
    )
