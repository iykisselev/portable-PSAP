#!/usr/bin/env python3
"""Append CADD, population AF, and Ensembl gene names to a cohort summary.

The Portable PSAP cohort summary deliberately contains only the candidate
signature and aggregate PSAP evidence.  This accessory command adds common
review fields without changing the certified workflow outputs.  It expects
the normal Portable PSAP ``candidate_signature`` spelling:
``CHROM:POS-POS_REF>ALT`` (with signatures separated by ``|``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path
import sys
from typing import TextIO

try:
    import pysam
except ImportError as exc:  # pragma: no cover - exercised by command users
    raise SystemExit(
        "This script requires pysam. Run it from the Portable PSAP validation "
        "environment or install pysam in the active environment."
    ) from exc


def open_text(path: Path, mode: str) -> TextIO:
    """Open a plain or gzip-compressed text path based on its magic bytes."""
    if "r" in mode and path.exists():
        with path.open("rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
    else:
        compressed = path.suffix == ".gz"
    return (
        gzip.open(path, mode, encoding="utf-8", newline="")
        if compressed
        else path.open(mode, encoding="utf-8", newline="")
    )


def canonical_chromosome(value: str) -> str:
    """Normalize only chromosome spelling for indexed resource lookup."""
    chromosome = value.strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    chromosome = chromosome.upper()
    return "MT" if chromosome == "M" else chromosome


def contig_lookup(contigs: tuple[str, ...] | list[str]) -> dict[str, str]:
    """Map canonical chromosome names to their unique resource spelling."""
    lookup: dict[str, str] = {}
    for contig in contigs:
        canonical = canonical_chromosome(contig)
        if canonical in lookup and lookup[canonical] != contig:
            raise ValueError(
                "Resource contains ambiguous contig spellings for "
                f"{canonical!r}: {lookup[canonical]!r}, {contig!r}"
            )
        lookup[canonical] = contig
    return lookup


def normalized_gene_id(value: str) -> str:
    """Discard the optional Ensembl version suffix before GTF lookup."""
    return value.strip().split(".", 1)[0]


def gtf_attributes(text: str) -> dict[str, str]:
    """Parse the simple key/value attribute field used by GTF files."""
    attributes: dict[str, str] = {}
    for raw_attribute in text.rstrip("\n").split(";"):
        fields = raw_attribute.strip().split(None, 1)
        if len(fields) == 2:
            attributes[fields[0]] = fields[1].strip().strip('"')
    return attributes


def read_gene_names(path: Path) -> dict[str, str]:
    """Return ``gene_id -> gene_name`` from gene records in a GTF."""
    names: dict[str, str] = {}
    with open_text(path, "rt") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}: GTF line {line_number} needs nine columns")
            if fields[2] != "gene":
                continue
            attributes = gtf_attributes(fields[8])
            gene_id = normalized_gene_id(attributes.get("gene_id", ""))
            gene_name = attributes.get("gene_name", "")
            if gene_id and gene_name:
                previous = names.setdefault(gene_id, gene_name)
                if previous != gene_name:
                    raise ValueError(
                        f"{path}: gene {gene_id!r} has conflicting names "
                        f"{previous!r} and {gene_name!r}"
                    )
    return names


def parse_variant_key(value: str) -> tuple[str, int, str, str]:
    """Parse one Portable PSAP ``CHROM:POS-POS_REF>ALT`` key."""
    try:
        chromosome, remainder = value.strip().split(":", 1)
        position_range, alleles = remainder.split("_", 1)
        position, end = position_range.split("-", 1)
        reference, alternate = alleles.split(">", 1)
        position_integer = int(position)
        if (
            position_integer < 1
            or int(end) < position_integer
            or not reference
            or not alternate
        ):
            raise ValueError
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid variant key {value!r}; expected CHROM:POS-POS_REF>ALT"
        ) from exc
    return (
        canonical_chromosome(chromosome),
        position_integer,
        reference.upper(),
        alternate.upper(),
    )


def format_value(value: object | None) -> str:
    """Render numeric lookup values while retaining a blank missing value."""
    if value is None:
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    # pysam decodes VCF Float INFO values through binary floating point.  Six
    # significant digits retain useful AF precision without exposing artifacts
    # such as ``0.0010000000474974513`` in a review table.
    return "" if not math.isfinite(number) else format(number, ".6g")


class AnnotationLookup:
    """Cached exact-allele lookups against CADD tables and a population VCF."""

    def __init__(
        self,
        cadd_snv: Path,
        cadd_indel: Path | None,
        population_vcf: Path,
        af_field: str,
    ) -> None:
        self.cadd_snv = pysam.TabixFile(str(cadd_snv))
        self.cadd_indel = pysam.TabixFile(str(cadd_indel)) if cadd_indel else None
        self.population = pysam.VariantFile(str(population_vcf))
        self.cadd_contigs = {
            "snv": contig_lookup(self.cadd_snv.contigs),
            "indel": contig_lookup(self.cadd_indel.contigs) if self.cadd_indel else {},
        }
        self.population_contigs = contig_lookup(list(self.population.header.contigs))
        self.af_field = af_field
        definition = self.population.header.info.get(af_field)
        if (
            definition is None
            or str(definition.number) not in {"A", "R", "1"}
            or definition.type not in {"Float", "Integer"}
        ):
            self.close()
            raise ValueError(f"INFO/{af_field} must declare numeric Number=A, R or 1")
        self.af_number = str(definition.number)
        self.cache: dict[tuple[str, int, str, str], tuple[str, str]] = {}

    def close(self) -> None:
        self.cadd_snv.close()
        if self.cadd_indel:
            self.cadd_indel.close()
        self.population.close()

    def values(self, key: tuple[str, int, str, str]) -> tuple[str, str]:
        """Return exact CADD PHRED and configured population AF for one allele."""
        if key not in self.cache:
            self.cache[key] = (self._cadd(key), self._af(key))
        return self.cache[key]

    def _cadd(self, key: tuple[str, int, str, str]) -> str:
        chromosome, position, reference, alternate = key
        kind = "snv" if len(reference) == len(alternate) == 1 else "indel"
        table = self.cadd_snv if kind == "snv" else self.cadd_indel
        resource_chromosome = self.cadd_contigs[kind].get(chromosome)
        if table is None or resource_chromosome is None:
            return ""
        try:
            rows = table.fetch(resource_chromosome, position - 1, position)
            for raw in rows:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                if (
                    canonical_chromosome(fields[0]) == chromosome
                    and int(fields[1]) == position
                    and fields[2].upper() == reference
                    and fields[3].upper() == alternate
                ):
                    return format_value(fields[5])
        except (ValueError, OSError) as exc:
            raise ValueError(
                f"Cannot query CADD resource {table.filename!r}: {exc}"
            ) from exc
        return ""

    def _af(self, key: tuple[str, int, str, str]) -> str:
        chromosome, position, reference, alternate = key
        resource_chromosome = self.population_contigs.get(chromosome)
        if resource_chromosome is None:
            return ""
        try:
            records = self.population.fetch(resource_chromosome, position - 1, position)
            for record in records:
                if record.pos != position or record.ref.upper() != reference:
                    continue
                alternates = tuple(value.upper() for value in record.alts or ())
                if alternate not in alternates or self.af_field not in record.info:
                    continue
                value = record.info[self.af_field]
                index = alternates.index(alternate)
                values = value if isinstance(value, tuple) else (value,)
                if self.af_number == "A" and len(values) == len(alternates):
                    return format_value(values[index])
                if self.af_number == "R" and len(values) == len(alternates) + 1:
                    return format_value(values[index + 1])
                if self.af_number == "1" and len(alternates) == len(values) == 1:
                    return format_value(values[0])
                raise ValueError(
                    f"Invalid INFO/{self.af_field} cardinality at {chromosome}:{position}"
                )
        except (ValueError, OSError) as exc:
            raise ValueError(
                f"Cannot query population resource {self.population.filename!r}: {exc}"
            ) from exc
        return ""


def require_distinct_output(summary: Path, output: Path) -> None:
    """Reject identical files, including hardlinks, before opening an output."""
    if summary.resolve() == output.resolve() or (
        output.exists() and summary.samefile(output)
    ):
        raise ValueError("Output cannot overwrite the input summary")


def annotate_summary(
    summary: Path,
    output: Path,
    gene_names: dict[str, str],
    lookup: AnnotationLookup,
    unit_column: str,
    signature_column: str,
) -> int:
    """Write a row-preserving, post-annotated cohort summary and return its size."""
    require_distinct_output(summary, output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(summary, "rt") as source, open_text(output, "wt") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{summary}: summary has no header")
        required = {unit_column, signature_column}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"{summary}: missing required column(s): {', '.join(sorted(missing))}"
            )
        additions = ("gene_name", "cadd_phred", "population_af")
        duplicate = set(additions).intersection(reader.fieldnames)
        if duplicate:
            raise ValueError(
                f"{summary}: annotation column(s) already exist: {', '.join(sorted(duplicate))}"
            )
        writer = csv.DictWriter(
            destination,
            fieldnames=[*reader.fieldnames, *additions],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        count = 0
        for row_number, row in enumerate(reader, start=2):
            signature = row[signature_column].strip()
            try:
                keys = [
                    parse_variant_key(value) for value in signature.split("|") if value
                ]
            except ValueError as exc:
                raise ValueError(f"{summary}: line {row_number}: {exc}") from exc
            if not keys:
                raise ValueError(
                    f"{summary}: line {row_number}: empty candidate signature"
                )
            cadd_values, af_values = zip(
                *(lookup.values(key) for key in keys), strict=True
            )
            row["gene_name"] = gene_names.get(normalized_gene_id(row[unit_column]), "")
            row["cadd_phred"] = ",".join(cadd_values)
            row["population_af"] = ",".join(af_values)
            writer.writerow(row)
            count += 1
    return count


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "summary", type=Path, help="Portable PSAP cohort_summary.tsv[.gz]"
    )
    parser.add_argument("output", type=Path, help="Annotated TSV or TSV.GZ output path")
    parser.add_argument(
        "--cadd-snv", required=True, type=Path, help="BGZF/tabix-indexed CADD SNV table"
    )
    parser.add_argument(
        "--cadd-indel", type=Path, help="BGZF/tabix-indexed CADD indel table"
    )
    parser.add_argument(
        "--population-vcf", required=True, type=Path, help="Indexed population VCF/BCF"
    )
    parser.add_argument(
        "--af-field",
        default="AF",
        help="Allele-frequency INFO field in the population VCF (default: AF)",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        type=Path,
        help="GTF or GTF.GZ containing Ensembl gene_id and gene_name",
    )
    parser.add_argument(
        "--unit-column",
        default="unit",
        help="Ensembl gene-ID column in the summary (default: unit)",
    )
    parser.add_argument(
        "--signature-column",
        default="candidate_signature",
        help="Variant-signature column (default: candidate_signature)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        require_distinct_output(args.summary, args.output)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    gene_names = read_gene_names(args.gtf)
    lookup = AnnotationLookup(
        args.cadd_snv, args.cadd_indel, args.population_vcf, args.af_field
    )
    try:
        rows = annotate_summary(
            args.summary,
            args.output,
            gene_names,
            lookup,
            args.unit_column,
            args.signature_column,
        )
    finally:
        lookup.close()
    print(f"Annotated {rows} cohort-summary rows: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
