"""Reference sequence and testing-unit membership construction."""

from __future__ import annotations

from bisect import bisect_right
from pathlib import Path

import polars as pl

from generation_partitions import GenomePartition, read_fasta_index
from genomic_contracts import canonical_chromosome

try:
    import pysam
except ImportError:
    pysam = None


class ReferenceReader:
    """Random-access FASTA reader using pysam or uncompressed FAI offsets."""

    def __init__(self, fasta: str | Path, fai: str | Path) -> None:
        self.fasta = Path(fasta)
        self.fai = Path(fai)
        self.contigs = read_fasta_index(self.fai)
        self.reader = pysam.FastaFile(str(self.fasta)) if pysam is not None else None
        self.layout: dict[str, tuple[int, int, int]] = {}
        if self.reader is None:
            with self.fai.open(encoding="utf-8") as handle:
                for raw in handle:
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        raise ValueError("Plain FASTA fallback needs five-column FAI records")
                    self.layout[canonical_chromosome(fields[0])] = (int(fields[2]), int(fields[3]), int(fields[4]))

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        """Fetch one bounded half-open interval with canonical chromosome names."""
        source_name, length = self.contigs[chromosome]
        if start < 0 or end > length or start > end:
            raise ValueError(f"Reference request is outside {chromosome!r}: {start}-{end}")
        if self.reader is not None:
            return self.reader.fetch(source_name, start, end).upper()

        offset, bases_per_line, bytes_per_line = self.layout[chromosome]
        sequence: list[bytes] = []
        position = start
        with self.fasta.open("rb") as handle:
            while position < end:
                within_line = position % bases_per_line
                take = min(end - position, bases_per_line - within_line)
                handle.seek(offset + (position // bases_per_line) * bytes_per_line + within_line)
                sequence.append(handle.read(take))
                position += take
        return b"".join(sequence).decode("ascii").upper()

    def close(self) -> None:
        """Close the indexed reader when pysam owns an open file handle."""
        if self.reader is not None:
            self.reader.close()


def membership_rows(
    partition: GenomePartition,
    mask: list[tuple[int, int]],
    intervals: list[tuple[str, int, int]],
    reference: ReferenceReader,
) -> pl.DataFrame:
    """Generate bounded candidate SNVs and distinct allele-to-unit memberships."""
    records: list[dict[str, object]] = []
    mask_starts = [start for start, _ in mask]
    for unit, start, end in intervals:
        mask_index = max(0, bisect_right(mask_starts, start) - 1)
        while mask_index < len(mask) and mask[mask_index][0] < end:
            mask_start, mask_end = mask[mask_index]
            overlap_start, overlap_end = max(start, mask_start), min(end, mask_end)
            mask_index += 1
            if overlap_start >= overlap_end:
                continue
            sequence = reference.fetch(partition.chromosome, overlap_start, overlap_end)
            for offset, base in enumerate(sequence, start=overlap_start):
                if base not in "ACGT":
                    continue
                for alternate in "ACGT":
                    if alternate != base:
                        records.append(
                            {
                                "testing_unit_id": unit,
                                "CHROM": partition.chromosome,
                                "POS": offset + 1,
                                "REF": base,
                                "ALT": alternate,
                            }
                        )
    return pl.DataFrame(
        records,
        schema={"testing_unit_id": pl.String, "CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String},
        strict=False,
    ).unique(maintain_order=True)
