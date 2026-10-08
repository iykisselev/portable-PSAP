#!/usr/bin/env python3
"""Validate PSAP workflow inputs before expensive rules run."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from bgzf import has_bgzf_eof, is_bgzf
from exclusions import data_rows, parse_exclusion_tsv_row
from genomic_contracts import canonical_chromosome, read_fai
from io_contracts import BedFormatError, is_missing_token, iter_bed, open_text
from sample_metadata import iter_sample_rows, normalize_sample

CANONICAL_CHROMOSOMES = tuple(str(chromosome) for chromosome in range(1, 23)) + (
    "X",
    "Y",
)
CANONICAL_CHROMOSOME_SET = set(CANONICAL_CHROMOSOMES)
CHROMOSOME_ORDER = {
    chromosome: position
    for position, chromosome in enumerate(CANONICAL_CHROMOSOMES)
}


@dataclass(frozen=True)
class IndexedResource:
    """A resource path and its optional index, kept together by contract."""

    path: str | Path
    index: str | Path | None = None


@dataclass(frozen=True)
class VcfHeaderFacts:
    """VCF header and contig facts shared by input validation checks."""

    samples: tuple[str, ...]
    contigs: frozenset[str]
    canonical_map: dict[str, str]


def genomic_resource_kind(path: str | Path) -> str | None:
    name = str(path).lower()
    for suffix in (".gz", ".bgz", ".bgzf"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for kind in ("vcf", "bed", "tsv", "txt"):
        if name.endswith(f".{kind}"):
            return "tsv" if kind == "txt" else kind
    return None


def require_file(
    path: str | Path,
    description: str,
    errors: list[str],
) -> bool:
    file_path = Path(path)
    if not file_path.is_file():
        errors.append(f"{description} does not exist or is not a file: {file_path}")
        return False

    if not os.access(file_path, os.R_OK):
        errors.append(f"{description} is not readable: {file_path}")
        return False

    return True


def validate_declared_index(
    index_path: str | Path | None,
    description: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    if is_missing_token(index_path):
        warnings.append(
            f"{description} has no configured index; the workflow will prepare "
            "an index under its output directory"
        )
        return
    require_file(index_path, f"{description} index", errors)


def open_text_maybe_gzip(path: str | Path):
    """Compatibility wrapper for validation's diagnostic decoding policy."""
    return open_text(path, errors="replace")


def has_data_row(
    path: str | Path,
    description: str,
    errors: list[str],
) -> bool:
    if not require_file(path, description, errors):
        return False

    try:
        with open_text_maybe_gzip(path) as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not any(value.strip() for value in row):
                    continue
                if row[0].lstrip().startswith("#"):
                    continue
                return True
    except OSError as exc:
        errors.append(f"Could not read {description}: {path}: {exc}")
        return False

    errors.append(f"No data rows found in {description}: {path}")
    return False


def validate_indexed_table(
    path: str | Path,
    description: str,
    errors: list[str],
    warnings: list[str],
    index_path: str | Path | None = None,
) -> None:
    if has_data_row(path, description, errors):
        bgzip = is_bgzf(path)
        if not bgzip:
            warnings.append(
                f"{description} is not BGZF-compressed; the workflow will prepare "
                f"a compressed copy under its output directory (source: {path})"
            )
        elif not has_bgzf_eof(path):
            warnings.append(
                f"{description} has no BGZF EOF marker; the workflow will prepare "
                f"a complete compressed copy under its output directory (source: {path})"
            )
        validate_declared_index(index_path, description, errors, warnings)


def validate_bed_regions(
    path: str | Path,
    description: str,
    errors: list[str],
    required_columns: int = 3,
    require_unit_id: bool = False,
) -> set[str]:
    if not require_file(path, description, errors):
        return set()

    data_rows = 0
    chromosomes: set[str] = set()
    def report(error: BedFormatError) -> None:
        nonlocal data_rows
        data_rows += 1
        errors.append(f"{description} {error}")

    try:
        for row in iter_bed(
            path, required_columns, errors="replace", on_error=report
        ):
            data_rows += 1
            if require_unit_id and is_missing_token(row.fields[3]):
                errors.append(
                    f"{description} line {row.line_number} has a missing unit ID "
                    "in column 4"
                )
            chromosomes.add(row.chromosome)
    except OSError as exc:
        errors.append(f"Could not read {description}: {path}: {exc}")
        return set()

    if data_rows == 0:
        errors.append(f"No data rows found in {description}: {path}")

    return chromosomes


def validate_exclusion_tsv(
    path: str | Path,
    description: str,
    errors: list[str],
) -> set[str]:
    """Validate 1-based interval TSV or CHROM/POS/REF/ALT variant TSV rows."""
    if not require_file(path, description, errors):
        return set()

    chromosomes: set[str] = set()
    data_count = 0
    try:
        for line_number, columns in data_rows(Path(path)):
            data_count += 1
            try:
                target = parse_exclusion_tsv_row(columns)
            except ValueError as exc:
                errors.append(f"{description} line {line_number} {exc}")
                continue
            chromosomes.add(target.chromosome)
    except OSError as exc:
        errors.append(f"Could not read {description}: {path}: {exc}")

    if data_count == 0:
        errors.append(f"No data rows found in {description}: {path}")
    return chromosomes


def validate_region_chromosomes(
    region_chromosomes: set[str],
    vcf_chromosomes: dict[str, str],
    description: str,
    warnings: list[str],
) -> None:
    if not region_chromosomes or not vcf_chromosomes:
        return

    vcf_canonical = set(vcf_chromosomes)
    missing = sorted(
        name
        for name in region_chromosomes
        if canonical_chromosome_name(name) not in vcf_canonical
    )
    if missing:
        warnings.append(
            f"{description} contains chromosome name(s) absent from the input VCF: "
            + ", ".join(missing)
        )


def validate_annotation_chromosome_names(
    resource_chromosomes: set[str],
    vcf_chromosomes: dict[str, str],
    description: str,
    warnings: list[str],
) -> None:
    if not resource_chromosomes or not vcf_chromosomes:
        return

    resource_names: dict[str, str] = {}
    for name in sorted(resource_chromosomes):
        canonical = canonical_chromosome_name(name)
        if canonical is not None and canonical not in resource_names:
            resource_names[canonical] = name

    mismatches = [
        (canonical, resource_names[canonical], vcf_chromosomes[canonical])
        for canonical in sorted(
            set(resource_names) & set(vcf_chromosomes),
            key=chromosome_sort_key,
        )
        if resource_names[canonical] != vcf_chromosomes[canonical]
    ]
    if mismatches:
        warnings.append(
            f"{description} chromosome names will be harmonized to the input VCF: "
            + ", ".join(
                f"{canonical}: resource={resource}, VCF={vcf}"
                for canonical, resource, vcf in mismatches
            )
        )


def validate_cadd_chromosome_names(
    resource_chromosomes: set[str],
    vcf_chromosomes: dict[str, str],
    description: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    """Require a canonical CADD/VCF contig overlap and report alias rewrites."""
    if not resource_chromosomes or not vcf_chromosomes:
        return

    resource_canonical = {
        canonical
        for name in resource_chromosomes
        if (canonical := canonical_chromosome_name(name)) is not None
    }
    shared = resource_canonical & set(vcf_chromosomes)
    if not shared:
        errors.append(
            f"{description} has no canonical chromosome overlap with the input "
            "VCF; check the resource assembly and chromosome names"
        )
        return

    validate_annotation_chromosome_names(
        resource_chromosomes,
        vcf_chromosomes,
        description,
        warnings,
    )


def validate_vcf_resource(
    path: str | Path,
    description: str,
    errors: list[str],
    warnings: list[str],
    index_path: str | Path | None = None,
) -> None:
    data_path = Path(path)
    name = str(data_path).lower()
    if not any(
        name.endswith(suffix)
        for suffix in (".vcf", ".vcf.gz", ".vcf.bgz", ".vcf.bgzf")
    ):
        errors.append(f"{description} must be a VCF file: {data_path}")
        return

    if not has_data_row(data_path, description, errors):
        return

    bgzip = is_bgzf(data_path)
    if not bgzip:
        warnings.append(
            f"{description} is not BGZF-compressed; the workflow will prepare a "
            f"compressed copy under its output directory (source: {data_path})"
        )
    elif not has_bgzf_eof(data_path):
        warnings.append(
            f"{description} has no BGZF EOF marker; the workflow will prepare a "
            f"complete compressed copy under its output directory (source: {data_path})"
        )

    validate_declared_index(index_path, description, errors, warnings)


def duplicate_values(values: list[str]) -> list[str]:
    return sorted(value for value, count in Counter(values).items() if count > 1)


def empty_ped_summary() -> dict:
    return {
        "manifest": [],
        "sample_ids": [],
        "all_sample_ids": [],
        "excluded_phenotype_sample_ids": [],
        "sex_counts": {},
        "phenotype_counts": {},
        "affected_count": 0,
        "male_count": 0,
    }


class _SampleMetadataAggregation:
    """Accumulate shared sample-metadata semantics for one input format."""

    def __init__(
        self,
        *,
        row_label: str,
        source_label: str,
        summary_label: str,
        duplicate_column: str,
        excluded_sample_label: str,
        sex_label: str,
        errors: list[str],
        warnings: list[str],
    ) -> None:
        self.row_label = row_label
        self.source_label = source_label
        self.summary_label = summary_label
        self.duplicate_column = duplicate_column
        self.excluded_sample_label = excluded_sample_label
        self.sex_label = sex_label
        self.errors = errors
        self.warnings = warnings
        self.all_sample_ids: list[str] = []
        self.sample_ids: list[str] = []
        self.excluded_phenotype_sample_ids: list[str] = []
        self.manifest: list[dict[str, str]] = []
        self.sex_values: list[str] = []
        self.phenotype_values: list[str] = []

    def add_row(
        self,
        line_number: int,
        family_id: str,
        sample_id: str,
        sex: str,
        phenotype: str,
    ) -> None:
        if not family_id:
            self.errors.append(
                f"{self.row_label} line {line_number} has an empty family ID"
            )

        if sample_id:
            self.all_sample_ids.append(sample_id)
            if "/" in sample_id or "\\" in sample_id:
                self.errors.append(
                    f"{self.row_label} line {line_number} sample ID contains a path "
                    f"separator: {sample_id!r}"
                )
            phenotype_label = {"1": "unaffected", "2": "affected"}.get(phenotype)
            if phenotype_label is None:
                self.excluded_phenotype_sample_ids.append(sample_id)
            else:
                self.sample_ids.append(sample_id)
                self.sex_values.append(sex)
                self.manifest.append(
                    {
                        "family_id": family_id or "NA",
                        "pid": sample_id,
                        "sex": sex if sex in {"1", "2"} else "NA",
                        "phenotype": phenotype_label,
                    }
                )
        else:
            self.errors.append(
                f"{self.row_label} line {line_number} has an empty sample ID"
            )

        self.phenotype_values.append(phenotype)

    def finalize(self, data_rows: int, path: str | Path) -> dict:
        if data_rows == 0:
            self.errors.append(f"No data rows found in {self.source_label}: {path}")

        duplicated = duplicate_values(self.all_sample_ids)
        if duplicated:
            self.errors.append(
                f"Duplicate sample IDs found in {self.duplicate_column}: "
                + ", ".join(duplicated)
            )

        phenotype_counts = Counter(self.phenotype_values)
        ignored_phenotypes = sorted(
            value for value in phenotype_counts if value not in {"1", "2"}
        )
        if self.excluded_phenotype_sample_ids and self.sample_ids:
            examples = ", ".join(self.excluded_phenotype_sample_ids[:10])
            remainder = len(self.excluded_phenotype_sample_ids) - min(
                len(self.excluded_phenotype_sample_ids),
                10,
            )
            suffix = f" (and {remainder} more)" if remainder else ""
            self.warnings.append(
                f"Excluded {len(self.excluded_phenotype_sample_ids)} "
                f"{self.excluded_sample_label} sample(s) with missing or unsupported "
                "phenotype values from all processing"
                + (
                    "; observed values: " + ", ".join(ignored_phenotypes)
                    if ignored_phenotypes
                    else ""
                )
                + f"; sample IDs: {examples}{suffix}"
            )

        sex_counts = Counter(self.sex_values)
        ignored_sex = sorted(value for value in sex_counts if value not in {"1", "2"})
        if ignored_sex:
            self.warnings.append(
                f"{self.sex_label} sex values other than 1 and 2 will not receive "
                "sex-specific chromosome-X recoding: "
                + ", ".join(ignored_sex)
            )

        if self.all_sample_ids and not self.sample_ids:
            self.errors.append(
                f"{self.summary_label} contains no samples with a known "
                "phenotype (1=unaffected or 2=affected)"
            )
        elif self.sample_ids:
            if phenotype_counts.get("2", 0) == 0:
                self.errors.append(
                    f"{self.summary_label} contains no affected individuals "
                    "(phenotype code 2)"
                )
            if phenotype_counts.get("1", 0) == 0:
                self.warnings.append(
                    f"{self.summary_label} contains no unaffected individuals "
                    "(phenotype code 1); the cohort report will use "
                    "validation=no_controls"
                )

        return {
            "manifest": self.manifest,
            "sample_ids": self.sample_ids,
            "all_sample_ids": self.all_sample_ids,
            "excluded_phenotype_sample_ids": self.excluded_phenotype_sample_ids,
            "sex_counts": dict(sorted(sex_counts.items())),
            "phenotype_counts": dict(sorted(phenotype_counts.items())),
            "affected_count": phenotype_counts.get("2", 0),
            "male_count": sex_counts.get("1", 0),
        }


def read_ped(
    path: str | Path,
    errors: list[str],
    warnings: list[str],
) -> dict:
    if not require_file(path, "PED file", errors):
        return empty_ped_summary()

    summary = _SampleMetadataAggregation(
        row_label="PED",
        source_label="PED file",
        summary_label="PED file",
        duplicate_column="PED column 2",
        excluded_sample_label="PED",
        sex_label="PED",
        errors=errors,
        warnings=warnings,
    )
    data_rows = 0

    try:
        def report(issue):
            nonlocal data_rows
            data_rows += 1
            errors.append(f"PED line {issue.line_number} {issue}")

        for row in iter_sample_rows(path, "ped", errors="replace", on_error=report):
            data_rows += 1
            sample = normalize_sample(row)
            summary.add_row(
                row.line_number,
                sample.family_id,
                sample.pid,
                sample.sex,
                sample.phenotype,
            )
    except OSError as exc:
        errors.append(f"Could not read PED file {path}: {exc}")

    return summary.finalize(data_rows, path)


def read_samples_tsv(
    path: str | Path,
    errors: list[str],
    warnings: list[str],
) -> dict:
    """Read the schema-validated TSV alternative to PED without changing coding."""
    if not require_file(path, "sample sheet", errors):
        return empty_ped_summary()

    summary = _SampleMetadataAggregation(
        row_label="Sample sheet",
        source_label="sample sheet",
        summary_label="Sample sheet",
        duplicate_column="sample sheet pid column",
        excluded_sample_label="sample-sheet",
        sex_label="Sample-sheet",
        errors=errors,
        warnings=warnings,
    )
    data_rows = 0
    try:
        header_error = False

        def report(issue):
            nonlocal header_error
            header_error = True
            errors.append(f"Sample sheet {issue}")

        for row in iter_sample_rows(path, "samples", on_error=report):
            data_rows += 1
            sample = normalize_sample(row)
            summary.add_row(
                row.line_number,
                sample.family_id,
                sample.pid,
                sample.sex,
                sample.phenotype,
            )
        if header_error:
            return empty_ped_summary()
    except OSError as exc:
        errors.append(f"Could not read sample sheet {path}: {exc}")

    return summary.finalize(data_rows, path)


def canonical_chromosome_name(name: str) -> str | None:
    canonical = canonical_chromosome(name)
    return canonical if canonical in CANONICAL_CHROMOSOME_SET else None


def chromosome_sort_key(chromosome: str) -> int:
    return CHROMOSOME_ORDER[chromosome]


def add_chromosome_name(
    chromosomes: dict[str, str],
    name: str,
    source: str,
    errors: list[str],
) -> None:
    canonical = canonical_chromosome_name(name)
    if canonical is None:
        return

    existing = chromosomes.get(canonical)
    if existing is None:
        chromosomes[canonical] = name
    elif existing != name:
        errors.append(
            f"{source} uses multiple names for canonical chromosome {canonical}: "
            f"{existing!r} and {name!r}"
        )


def parse_contig_id(header_line: str) -> str | None:
    prefix = "##contig=<ID="
    if not header_line.startswith(prefix):
        return None

    contig_id = header_line[len(prefix) :].split(",", 1)[0].split(">", 1)[0]
    return contig_id.strip() or None


def read_vcf_record_contigs(
    path: str | Path,
    errors: list[str],
) -> set[str]:
    contigs: set[str] = set()

    try:
        with open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                name = line.split("\t", 1)[0].strip()
                if name:
                    contigs.add(name)
    except OSError as exc:
        errors.append(f"Could not read chromosome names in VCF {path}: {exc}")

    return contigs


def read_indexed_vcf_contigs(path: str | Path, bcftools: str) -> set[str]:
    """Read indexed VCF contigs without decompressing records, when possible."""
    try:
        result = subprocess.run(
            [bcftools, "index", "--stats", str(path)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError:
        return set()
    if result.returncode != 0:
        return set()
    return {
        line.split("\t", 1)[0].strip()
        for line in result.stdout.splitlines()
        if line.strip()
    }


def read_vcf_header_facts(path: str | Path, errors: list[str]) -> VcfHeaderFacts:
    """Read VCF samples and contigs, using an index before a record scan."""
    if not require_file(path, "Input VCF", errors):
        return VcfHeaderFacts((), frozenset(), {})

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        errors.append("bcftools is not available")
        return VcfHeaderFacts((), frozenset(), {})

    try:
        result = subprocess.run(
            [bcftools, "view", "-h", str(path)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError as exc:
        errors.append(f"Could not run bcftools to read VCF {path}: {exc}")
        return VcfHeaderFacts((), frozenset(), {})

    if result.returncode != 0:
        errors.append(
            f"bcftools could not read the VCF header: {path}\n"
            f"bcftools stderr: {result.stderr.strip()}"
        )
        return VcfHeaderFacts((), frozenset(), {})

    samples: tuple[str, ...] = ()
    contigs: set[str] = set()
    for line in result.stdout.splitlines():
        contig_id = parse_contig_id(line)
        if contig_id is not None:
            contigs.add(contig_id)
        elif line.startswith("#CHROM"):
            samples = tuple(sample.strip() for sample in line.split("\t")[9:])

    if not contigs:
        contigs = read_indexed_vcf_contigs(path, bcftools)
    if not contigs:
        contigs = read_vcf_record_contigs(path, errors)

    chromosomes: dict[str, str] = {}
    for contig in sorted(contigs):
        add_chromosome_name(chromosomes, contig, "Input VCF", errors)
    if not chromosomes:
        errors.append("No canonical chromosomes (1-22, X, Y) were found in the VCF")

    return VcfHeaderFacts(samples, frozenset(contigs), chromosomes)


def read_vcf_record_chromosomes(
    path: str | Path,
    errors: list[str],
) -> dict[str, str]:
    """Compatibility wrapper for callers that need only canonical contigs."""
    chromosomes: dict[str, str] = {}
    for contig in sorted(read_vcf_record_contigs(path, errors)):
        add_chromosome_name(chromosomes, contig, "Input VCF", errors)
    return chromosomes


def read_vcf(
    path: str | Path,
    errors: list[str],
) -> tuple[list[str], dict[str, str]]:
    facts = read_vcf_header_facts(path, errors)
    return list(facts.samples), facts.canonical_map


def read_all_vcf_contigs(path: str | Path, bcftools: str) -> set[str]:
    result = subprocess.run(
        [bcftools, "view", "--header-only", str(path)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        return set()
    contigs = {
        contig
        for line in result.stdout.splitlines()
        if (contig := parse_contig_id(line)) is not None
    }
    if contigs:
        return contigs

    with open_text_maybe_gzip(path) as handle:
        return {
            line.split("\t", 1)[0]
            for line in handle
            if line and not line.startswith("#")
        }


def read_fasta_contigs(index_path: str | Path) -> dict[str, str]:
    index_path = Path(index_path)
    if not index_path.is_file():
        return {}
    return {chromosome: record.name for chromosome, record in read_fai(index_path).items()}


def validate_normalized_split_vcf(
    vcf_path: str | Path,
    fasta_path: str | Path,
    fasta_index_path: str | Path,
    errors: list[str],
    warnings: list[str],
    vcf_facts: VcfHeaderFacts | None = None,
) -> None:
    if not Path(vcf_path).is_file() or not Path(fasta_path).is_file():
        return

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        return

    vcf_contigs = (
        set(vcf_facts.contigs)
        if vcf_facts is not None
        else read_all_vcf_contigs(vcf_path, bcftools)
    )
    try:
        fasta_contigs = read_fasta_contigs(fasta_index_path)
    except ValueError as exc:
        errors.append(f"Reference FASTA index is invalid: {exc}")
        return
    vcf_by_canonical = {
        canonical_chromosome(contig): contig for contig in vcf_contigs
    }
    shared_contigs = sorted(set(vcf_by_canonical) & set(fasta_contigs))
    unsupported_contigs = sorted(
        contig for canonical, contig in vcf_by_canonical.items()
        if canonical not in fasta_contigs
    )
    if unsupported_contigs:
        examples = ", ".join(unsupported_contigs[:10])
        remainder = len(unsupported_contigs) - min(len(unsupported_contigs), 10)
        suffix = f" (and {remainder} more)" if remainder else ""
        warnings.append(
            "Input VCF contigs absent from the reference FASTA will be skipped "
            f"during normalization validation: {examples}{suffix}"
        )

    if not shared_contigs:
        errors.append(
            "Input VCF and reference FASTA have no canonically compatible contigs; "
            "check the assembly and chromosome naming convention"
        )


def validate_vcf_samples(
    vcf_samples: list[str],
    metadata_samples: list[str],
    errors: list[str],
    metadata_label: str = "PED",
) -> None:
    if not vcf_samples:
        errors.append("No sample IDs were found in the input VCF")
        return

    duplicated = duplicate_values(vcf_samples)
    if duplicated:
        errors.append(
            "Duplicate sample IDs found in VCF header: " + ", ".join(duplicated)
        )

    missing = sorted(set(metadata_samples) - set(vcf_samples))
    if missing:
        errors.append(
            f"{metadata_label} sample IDs missing from the VCF header: "
            + ", ".join(missing)
        )


def validate_hemizygote_chromosomes(
    hemizygote: bool,
    vcf_chromosomes: dict[str, str],
    errors: list[str],
) -> None:
    if not hemizygote or not vcf_chromosomes:
        return

    missing = [chromosome for chromosome in ("X",) if chromosome not in vcf_chromosomes]
    if missing:
        errors.append(
            "HEM chromosome-X model requires chromosome(s) in the input VCF: "
            + ", ".join(missing)
        )


def validate_lookup_files(
    lookup_paths: dict[str, str],
    settings: dict[str, bool],
    errors: list[str],
) -> list[str]:
    required = ["het", "hom"]
    if settings["hemizygote"]:
        required.append("hem_chrX")

    checked_paths: list[str] = []
    for name in required:
        path = lookup_paths.get(name, "").strip()
        if is_missing_token(path):
            errors.append(
                f"PSAP lookup path for {name!r} is required when its model is enabled"
            )
            continue
        has_data_row(path, f"PSAP {name} lookup file", errors)
        checked_paths.append(path)

    return checked_paths


def required_lookup_names(settings: dict[str, bool]) -> list[str]:
    names = ["het", "hom"]
    if settings["hemizygote"]:
        names.append("hem_chrX")
    return names


def read_testing_unit_ids(
    path: str | Path,
    errors: list[str],
) -> set[str]:
    if not require_file(path, "Testing-unit BED file", errors):
        return set()

    identifiers: set[str] = set()
    try:
        with open_text_maybe_gzip(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if (
                    not stripped
                    or stripped.startswith("#")
                    or stripped.startswith("track ")
                    or stripped.startswith("browser ")
                ):
                    continue
                columns = stripped.split("\t")
                if len(columns) < 4:
                    continue
                identifier = columns[3].strip()
                if not is_missing_token(identifier):
                    identifiers.add(identifier)
    except OSError as exc:
        errors.append(f"Could not read testing-unit IDs from {path}: {exc}")
    return identifiers


def read_lookup_unit_ids(
    path: str | Path,
    description: str,
    errors: list[str],
) -> set[str]:
    if not Path(path).is_file():
        return set()

    identifiers: set[str] = set()
    try:
        with open_text_maybe_gzip(path) as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                    continue
                identifier = row[0].strip()
                if identifier.lower() == "testing_unit_id":
                    continue
                identifiers.add(identifier)
    except OSError as exc:
        errors.append(f"Could not read testing-unit IDs from {description}: {path}: {exc}")
    return identifiers


def validate_lookup_unit_overlap(
    testing_unit_bed: str | Path,
    lookup_paths: dict[str, str],
    settings: dict[str, bool],
    errors: list[str],
) -> dict[str, dict[str, int]]:
    """Reject lookup sets that cannot possibly score the configured units."""
    testing_unit_ids = read_testing_unit_ids(testing_unit_bed, errors)
    summary: dict[str, dict[str, int]] = {}

    for name in required_lookup_names(settings):
        path = lookup_paths.get(name, "").strip()
        if is_missing_token(path):
            continue
        lookup_ids = read_lookup_unit_ids(
            path,
            f"PSAP {name} lookup file",
            errors,
        )
        overlap = testing_unit_ids & lookup_ids
        summary[name] = {
            "testing_unit_id_count": len(testing_unit_ids),
            "lookup_id_count": len(lookup_ids),
            "overlap_count": len(overlap),
        }
        if testing_unit_ids and lookup_ids and not overlap:
            errors.append(
                f"PSAP {name} lookup file has zero testing-unit ID overlap with "
                f"the configured testing-unit BED: {path}"
            )

    return summary


def validate_reference_tables(
    *,
    calibration_regions: IndexedResource | None,
    exclusions: Sequence[IndexedResource],
    testing_units: IndexedResource,
    par_regions: Sequence[IndexedResource],
    unit: str,
    errors: list[str],
    warnings: list[str],
) -> dict[str, set[str]]:
    if calibration_regions is None:
        warnings.append(
            "No generated calibration mask supplied; only configured exclusions will filter regions"
        )
    else:
        calibration_regions_kind = genomic_resource_kind(calibration_regions.path)
        if calibration_regions_kind not in {"bed", "tsv"}:
            errors.append(
                "Generated calibration regions must be a BED or TSV file; received "
                f"{calibration_regions.path!r}"
            )
        else:
            validate_indexed_table(
                calibration_regions.path,
                "Calibration regions file",
                errors,
                warnings,
                calibration_regions.index,
            )

    unit_description = f"Configured {unit} testing-unit BED file"
    validate_indexed_table(
        testing_units.path,
        unit_description,
        errors,
        warnings,
        testing_units.index,
    )
    unit_chromosomes = validate_bed_regions(
        testing_units.path,
        unit_description,
        errors,
        required_columns=4,
        require_unit_id=True,
    )

    configured_exclusions = [*exclusions, *par_regions]
    if not configured_exclusions:
        warnings.append(
            "Neither reference.variants_exclude nor reference.par_regions was "
            "provided; exclusion filtering will be skipped"
        )

    exclusion_chromosomes: set[str] = set()
    for position, resource in enumerate(configured_exclusions, start=1):
        path = resource.path
        description = f"Exclusion file {position}"
        kind = genomic_resource_kind(path)
        if kind is None:
            errors.append(
                f"{description} must use a .bed, .tsv, .txt, or .vcf extension: {path}"
            )
            continue
        if kind == "vcf":
            validate_vcf_resource(path, description, errors, warnings, resource.index)
            continue

        validate_indexed_table(path, description, errors, warnings, resource.index)
        if kind == "bed":
            exclusion_chromosomes.update(
                validate_bed_regions(path, description, errors)
            )
        else:
            exclusion_chromosomes.update(
                validate_exclusion_tsv(path, description, errors)
            )

    return {
        "testing_unit": unit_chromosomes,
        "exclusions": exclusion_chromosomes,
    }


def validate_fasta(
    path: str | Path,
    index_path: str | Path,
    gzi_path: str | Path | None,
    errors: list[str],
) -> None:
    if require_file(path, "Reference FASTA", errors):
        require_file(index_path, "Reference FASTA index", errors)
        if str(path).endswith(".gz"):
            if is_missing_token(gzi_path):
                errors.append("BGZF FASTA index binding is required for a compressed reference FASTA")
            else:
                require_file(gzi_path, "BGZF FASTA index", errors)


def validate_cadd_resources(
    cadd_snv: str | Path,
    cadd_indel: str | Path | None,
    include_indels: bool,
    errors: list[str],
    warnings: list[str],
    cadd_snv_index: str | Path | None = None,
    cadd_indel_index: str | Path | None = None,
) -> dict[str, set[str]]:
    resources = [("whole_genome_SNVs.tsv.gz", Path(cadd_snv), cadd_snv_index)]
    if include_indels:
        if is_missing_token(cadd_indel):
            errors.append("CADD indel resource is required when include_indels is enabled")
            return {}
        resources.append(("InDels.tsv.gz", Path(cadd_indel), cadd_indel_index))

    chromosomes: dict[str, set[str]] = {}
    for filename, path, index_path in resources:
        description = f"CADD resource {filename}"
        validate_indexed_table(path, description, errors, warnings, index_path)
        validate_cadd_columns(path, description, errors)
        chromosomes[filename] = read_indexed_table_chromosomes(
            path,
            description,
            errors,
            index_path,
        )

    return chromosomes


def read_indexed_table_chromosomes(
    path: str | Path,
    description: str,
    errors: list[str],
    index_path: str | Path | None = None,
    sample_rows: int = 10000,
) -> set[str]:
    """Read contigs from the bound index when tabix can consume it, else sample."""
    data_path = Path(path)
    if not data_path.is_file():
        return set()

    tabix = shutil.which("tabix")
    configured_index = Path(index_path) if index_path else None
    tabix_compatible_index = configured_index in (
        Path(f"{data_path}.tbi"),
        Path(f"{data_path}.csi"),
    )
    if tabix is not None and tabix_compatible_index:
        try:
            result = subprocess.run(
                [tabix, "--list-chroms", str(data_path)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
        except OSError as exc:
            errors.append(f"Could not inspect chromosome names in {description}: {exc}")
            return set()
        if result.returncode == 0:
            names = {line.strip() for line in result.stdout.splitlines() if line.strip()}
            if names:
                return names

    names: set[str] = set()
    inspected = 0
    try:
        with open_text_maybe_gzip(data_path) as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                names.add(stripped.split("\t", 1)[0].strip())
                inspected += 1
                if inspected >= sample_rows:
                    break
    except OSError as exc:
        errors.append(f"Could not inspect chromosome names in {description}: {exc}")
    return {name for name in names if name}


def validate_cadd_columns(
    path: str | Path,
    description: str,
    errors: list[str],
) -> None:
    if not Path(path).is_file():
        return

    try:
        with open_text_maybe_gzip(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                columns = stripped.split("\t")
                if len(columns) < 6:
                    errors.append(
                        f"{description} line {line_number} must contain CADD "
                        "columns Chrom, Pos, Ref, Alt, RawScore, and PHRED"
                    )
                    return

                try:
                    position = int(columns[1])
                    float(columns[5])
                except ValueError:
                    errors.append(
                        f"{description} line {line_number} has a non-numeric "
                        "position or PHRED score"
                    )
                    return

                if position < 1 or not columns[0] or not columns[2] or not columns[3]:
                    errors.append(
                        f"{description} line {line_number} has invalid variant "
                        "coordinates or alleles"
                    )
                return
    except OSError as exc:
        errors.append(f"Could not inspect {description}: {path}: {exc}")
