#!/usr/bin/env python3
"""Validate PSAP workflow inputs before expensive rules run."""

from __future__ import annotations

import csv
import gzip
import os
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from bgzf import has_bgzf_eof, is_bgzf
from genomic_contracts import canonical_chromosome

VALID_UNITS = {"gene", "cadd_region", "coding_cadd_region"}
BOOLEAN_VALUES = {"true": True, "false": False}
MISSING_CONFIG_VALUES = {"", "na", "none", "null", ".", "-"}
CANONICAL_CHROMOSOMES = tuple(str(chromosome) for chromosome in range(1, 23)) + (
    "X",
    "Y",
)
CANONICAL_CHROMOSOME_SET = set(CANONICAL_CHROMOSOMES)
CHROMOSOME_ORDER = {
    chromosome: position
    for position, chromosome in enumerate(CANONICAL_CHROMOSOMES)
}


def optional_config_paths(value: str) -> list[str]:
    normalized = value.strip()
    if normalized.lower() in MISSING_CONFIG_VALUES:
        return []
    return [
        item.strip()
        for item in normalized.split(",")
        if item.strip().lower() not in MISSING_CONFIG_VALUES
    ]


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


def parse_bool(value: str, name: str, errors: list[str]) -> bool:
    normalized = value.strip().lower()
    if normalized in BOOLEAN_VALUES:
        return BOOLEAN_VALUES[normalized]

    errors.append(f"{name} must be True or False; received {value!r}")
    return False


def parse_model_settings(
    args: object,
    errors: list[str],
) -> dict[str, bool]:
    if args.unit not in VALID_UNITS:
        errors.append(
            "unit must be one of "
            + ", ".join(sorted(VALID_UNITS))
            + f"; received {args.unit!r}"
        )

    if not args.assembly.strip():
        errors.append("assembly cannot be empty")

    validate_positive_number(args.score_max, "score_max", errors)
    settings = {
        "hemizygote": parse_bool(
            args.hemizygote_model,
            "hemizygote_model",
            errors,
        ),
        "include_indels": parse_bool(args.include_indels, "include_indels", errors),
    }

    return settings


def validate_positive_number(value: str, name: str, errors: list[str]) -> None:
    try:
        number = float(value)
    except ValueError:
        errors.append(f"{name} must be numeric; received {value!r}")
        return

    if number <= 0:
        errors.append(f"{name} must be positive; received {value!r}")


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
    if index_path is None or str(index_path).strip().lower() in MISSING_CONFIG_VALUES:
        warnings.append(
            f"{description} has no configured index; the workflow will prepare "
            "an index under its output directory"
        )
        return
    require_file(index_path, f"{description} index", errors)


def open_text_maybe_gzip(path: str | Path):
    file_path = Path(path)
    try:
        with file_path.open("rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
    except OSError:
        compressed = file_path.suffix in {".gz", ".bgz", ".bgzf"}
    if compressed:
        return gzip.open(file_path, "rt", encoding="utf-8", errors="replace")

    return file_path.open("rt", encoding="utf-8", errors="replace")


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

                data_rows += 1
                columns = stripped.split("\t")
                if len(columns) < required_columns:
                    errors.append(
                        f"{description} line {line_number} must contain at least "
                        f"{required_columns} tab-delimited BED columns"
                    )
                    continue

                if require_unit_id and columns[3].strip() in MISSING_CONFIG_VALUES:
                    errors.append(
                        f"{description} line {line_number} has a missing unit ID "
                        "in column 4"
                    )

                chromosome = columns[0].strip()
                if not chromosome:
                    errors.append(
                        f"{description} line {line_number} has an empty chromosome"
                    )
                else:
                    chromosomes.add(chromosome)

                try:
                    start = int(columns[1])
                    end = int(columns[2])
                except ValueError:
                    errors.append(
                        f"{description} line {line_number} has non-integer "
                        "start or end coordinates"
                    )
                    continue

                if start < 0 or end <= start:
                    errors.append(
                        f"{description} line {line_number} must satisfy "
                        f"0 <= start < end; received {start}, {end}"
                    )
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
        with open_text_maybe_gzip(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                data_count += 1
                columns = stripped.split("\t")
                if len(columns) < 3:
                    errors.append(
                        f"{description} line {line_number} must have at least 3 columns"
                    )
                    continue
                chromosome = columns[0].strip()
                if not chromosome:
                    errors.append(
                        f"{description} line {line_number} has an empty chromosome"
                    )
                else:
                    chromosomes.add(chromosome)
                try:
                    start = int(columns[1])
                except ValueError:
                    errors.append(
                        f"{description} line {line_number} has a non-integer position"
                    )
                    continue

                try:
                    end = int(columns[2])
                except ValueError:
                    if len(columns) < 4 or start < 1 or not columns[2] or not columns[3]:
                        errors.append(
                            f"{description} line {line_number} must contain either "
                            "CHROM,START,END or CHROM,POS,REF,ALT"
                        )
                else:
                    if start < 1 or end < start:
                        errors.append(
                            f"{description} line {line_number} interval must satisfy "
                            f"1 <= start <= end; received {start}, {end}"
                        )
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


def read_ped(
    path: str | Path,
    errors: list[str],
    warnings: list[str],
) -> dict:
    if not require_file(path, "PED file", errors):
        return empty_ped_summary()

    all_sample_ids: list[str] = []
    sample_ids: list[str] = []
    excluded_phenotype_sample_ids: list[str] = []
    manifest: list[dict[str, str]] = []
    sex_values: list[str] = []
    phenotype_values: list[str] = []
    data_rows = 0

    try:
        with Path(path).open(
            "rt",
            encoding="utf-8",
            errors="replace",
            newline="",
        ) as handle:
            for line_number, row in enumerate(
                csv.reader(handle, delimiter="\t"),
                start=1,
            ):
                if not row or all(not value.strip() for value in row):
                    continue

                data_rows += 1
                if len(row) < 6:
                    errors.append(
                        f"PED line {line_number} has {len(row)} columns; "
                        "at least 6 are required"
                    )
                    continue

                family_id = row[0].strip()
                sample_id = row[1].strip()
                sex = row[4].strip()
                phenotype = row[5].strip()

                if not family_id:
                    errors.append(f"PED line {line_number} has an empty family ID")

                if sample_id:
                    all_sample_ids.append(sample_id)
                    if "/" in sample_id or "\\" in sample_id:
                        errors.append(
                            f"PED line {line_number} sample ID contains a path "
                            f"separator: {sample_id!r}"
                        )
                    phenotype_label = {
                        "1": "unaffected",
                        "2": "affected",
                    }.get(phenotype)
                    if phenotype_label is None:
                        excluded_phenotype_sample_ids.append(sample_id)
                    else:
                        sample_ids.append(sample_id)
                        sex_values.append(sex)
                        manifest.append(
                            {
                                "family_id": family_id or "NA",
                                "pid": sample_id,
                                "sex": sex if sex in {"1", "2"} else "NA",
                                "phenotype": phenotype_label,
                            }
                        )
                else:
                    errors.append(f"PED line {line_number} has an empty sample ID")

                phenotype_values.append(phenotype)
    except OSError as exc:
        errors.append(f"Could not read PED file {path}: {exc}")

    if data_rows == 0:
        errors.append(f"No data rows found in PED file: {path}")

    duplicated = duplicate_values(all_sample_ids)
    if duplicated:
        errors.append(
            "Duplicate sample IDs found in PED column 2: " + ", ".join(duplicated)
        )

    phenotype_counts = Counter(phenotype_values)
    ignored_phenotypes = sorted(
        value for value in phenotype_counts if value not in {"1", "2"}
    )
    if excluded_phenotype_sample_ids and sample_ids:
        examples = ", ".join(excluded_phenotype_sample_ids[:10])
        remainder = len(excluded_phenotype_sample_ids) - min(
            len(excluded_phenotype_sample_ids),
            10,
        )
        suffix = f" (and {remainder} more)" if remainder else ""
        warnings.append(
            f"Excluded {len(excluded_phenotype_sample_ids)} PED sample(s) with "
            "missing or unsupported phenotype values from all processing"
            + (
                "; observed values: " + ", ".join(ignored_phenotypes)
                if ignored_phenotypes
                else ""
            )
            + f"; sample IDs: {examples}{suffix}"
        )

    sex_counts = Counter(sex_values)
    ignored_sex = sorted(value for value in sex_counts if value not in {"1", "2"})
    if ignored_sex:
        warnings.append(
            "PED sex values other than 1 and 2 will not receive sex-specific "
            "chromosome-X recoding: "
            + ", ".join(ignored_sex)
        )

    if all_sample_ids and not sample_ids:
        errors.append(
            "PED file contains no samples with a known phenotype "
            "(1=unaffected or 2=affected)"
        )
    elif sample_ids:
        if phenotype_counts.get("2", 0) == 0:
            errors.append("PED file contains no affected individuals (phenotype code 2)")

        if phenotype_counts.get("1", 0) == 0:
            warnings.append(
                "PED file contains no unaffected individuals (phenotype code 1); "
                "the cohort report will use validation=no_controls"
            )

    return {
        "manifest": manifest,
        "sample_ids": sample_ids,
        "all_sample_ids": all_sample_ids,
        "excluded_phenotype_sample_ids": excluded_phenotype_sample_ids,
        "sex_counts": dict(sorted(sex_counts.items())),
        "phenotype_counts": dict(sorted(phenotype_counts.items())),
        "affected_count": phenotype_counts.get("2", 0),
        "male_count": sex_counts.get("1", 0),
    }


def read_samples_tsv(
    path: str | Path,
    errors: list[str],
    warnings: list[str],
) -> dict:
    """Read the schema-validated TSV alternative to PED without changing coding."""
    if not require_file(path, "sample sheet", errors):
        return empty_ped_summary()

    all_sample_ids: list[str] = []
    sample_ids: list[str] = []
    excluded_phenotype_sample_ids: list[str] = []
    manifest: list[dict[str, str]] = []
    sex_values: list[str] = []
    phenotype_values: list[str] = []
    data_rows = 0
    required_columns = ["family_id", "pid", "sex", "phenotype"]

    try:
        with Path(path).open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers = reader.fieldnames or []
            missing = [name for name in required_columns if name not in headers]
            if missing:
                errors.append(
                    "Sample sheet is missing required column(s): " + ", ".join(missing)
                )
                return empty_ped_summary()

            for line_number, row in enumerate(reader, start=2):
                if not row or not any((value or "").strip() for value in row.values()):
                    continue
                data_rows += 1
                family_id = (row.get("family_id") or "").strip()
                sample_id = (row.get("pid") or "").strip()
                sex = (row.get("sex") or "").strip()
                phenotype = (row.get("phenotype") or "").strip()

                if not family_id:
                    errors.append(f"Sample sheet line {line_number} has an empty family ID")
                if not sample_id:
                    errors.append(f"Sample sheet line {line_number} has an empty sample ID")
                else:
                    all_sample_ids.append(sample_id)
                    if "/" in sample_id or "\\" in sample_id:
                        errors.append(
                            f"Sample sheet line {line_number} sample ID contains a path "
                            f"separator: {sample_id!r}"
                        )
                    phenotype_label = {"1": "unaffected", "2": "affected"}.get(phenotype)
                    if phenotype_label is None:
                        excluded_phenotype_sample_ids.append(sample_id)
                    else:
                        sample_ids.append(sample_id)
                        sex_values.append(sex)
                        manifest.append(
                            {
                                "family_id": family_id or "NA",
                                "pid": sample_id,
                                "sex": sex if sex in {"1", "2"} else "NA",
                                "phenotype": phenotype_label,
                            }
                        )
                phenotype_values.append(phenotype)
    except OSError as exc:
        errors.append(f"Could not read sample sheet {path}: {exc}")

    if data_rows == 0:
        errors.append(f"No data rows found in sample sheet: {path}")

    duplicated = duplicate_values(all_sample_ids)
    if duplicated:
        errors.append("Duplicate sample IDs found in sample sheet pid column: " + ", ".join(duplicated))

    phenotype_counts = Counter(phenotype_values)
    ignored_phenotypes = sorted(
        value for value in phenotype_counts if value not in {"1", "2"}
    )
    if excluded_phenotype_sample_ids and sample_ids:
        examples = ", ".join(excluded_phenotype_sample_ids[:10])
        remainder = len(excluded_phenotype_sample_ids) - min(
            len(excluded_phenotype_sample_ids), 10
        )
        suffix = f" (and {remainder} more)" if remainder else ""
        warnings.append(
            f"Excluded {len(excluded_phenotype_sample_ids)} sample-sheet sample(s) "
            "with missing or unsupported phenotype values from all processing"
            + (
                "; observed values: " + ", ".join(ignored_phenotypes)
                if ignored_phenotypes
                else ""
            )
            + f"; sample IDs: {examples}{suffix}"
        )

    sex_counts = Counter(sex_values)
    ignored_sex = sorted(value for value in sex_counts if value not in {"1", "2"})
    if ignored_sex:
        warnings.append(
            "Sample-sheet sex values other than 1 and 2 will not receive sex-specific "
            "chromosome-X recoding: " + ", ".join(ignored_sex)
        )

    if all_sample_ids and not sample_ids:
        errors.append(
            "Sample sheet contains no samples with a known phenotype "
            "(1=unaffected or 2=affected)"
        )
    elif sample_ids:
        if phenotype_counts.get("2", 0) == 0:
            errors.append("Sample sheet contains no affected individuals (phenotype code 2)")
        if phenotype_counts.get("1", 0) == 0:
            warnings.append(
                "Sample sheet contains no unaffected individuals (phenotype code 1); "
                "the cohort report will use validation=no_controls"
            )

    return {
        "manifest": manifest,
        "sample_ids": sample_ids,
        "all_sample_ids": all_sample_ids,
        "excluded_phenotype_sample_ids": excluded_phenotype_sample_ids,
        "sex_counts": dict(sorted(sex_counts.items())),
        "phenotype_counts": dict(sorted(phenotype_counts.items())),
        "affected_count": phenotype_counts.get("2", 0),
        "male_count": sex_counts.get("1", 0),
    }


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


def read_vcf_record_chromosomes(
    path: str | Path,
    errors: list[str],
) -> dict[str, str]:
    chromosomes: dict[str, str] = {}
    seen_names: set[str] = set()

    try:
        with open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                name = line.split("\t", 1)[0].strip()
                if name and name not in seen_names:
                    seen_names.add(name)
                    add_chromosome_name(chromosomes, name, "Input VCF", errors)
    except OSError as exc:
        errors.append(f"Could not read chromosome names in VCF {path}: {exc}")

    return chromosomes


def read_vcf(
    path: str | Path,
    errors: list[str],
) -> tuple[list[str], dict[str, str]]:
    if not require_file(path, "Input VCF", errors):
        return [], {}

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        errors.append("bcftools is not available")
        return [], {}

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
        return [], {}

    if result.returncode != 0:
        errors.append(
            f"bcftools could not read the VCF header: {path}\n"
            f"bcftools stderr: {result.stderr.strip()}"
        )
        return [], {}

    samples: list[str] = []
    chromosomes: dict[str, str] = {}
    for line in result.stdout.splitlines():
        contig_id = parse_contig_id(line)
        if contig_id is not None:
            add_chromosome_name(chromosomes, contig_id, "Input VCF", errors)
        elif line.startswith("#CHROM"):
            samples = [sample.strip() for sample in line.split("\t")[9:]]

    if not chromosomes:
        chromosomes = read_vcf_record_chromosomes(path, errors)

    if not chromosomes:
        errors.append("No canonical chromosomes (1-22, X, Y) were found in the VCF")

    return samples, chromosomes


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


def read_fasta_contigs(index_path: str | Path) -> set[str]:
    index_path = Path(index_path)
    if not index_path.is_file():
        return set()
    with index_path.open("rt", encoding="utf-8", errors="replace") as handle:
        return {
            line.split("\t", 1)[0].strip()
            for line in handle
            if line.strip()
        }


def validate_normalized_split_vcf(
    vcf_path: str | Path,
    fasta_path: str | Path,
    fasta_index_path: str | Path,
    errors: list[str],
    warnings: list[str],
) -> None:
    if not Path(vcf_path).is_file() or not Path(fasta_path).is_file():
        return

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        return

    vcf_contigs = read_all_vcf_contigs(vcf_path, bcftools)
    fasta_contigs = read_fasta_contigs(fasta_index_path)
    shared_contigs = sorted(vcf_contigs & fasta_contigs)
    unsupported_contigs = sorted(vcf_contigs - fasta_contigs)
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
            "Input VCF and reference FASTA have no identically named contigs; "
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
        if path.lower() in MISSING_CONFIG_VALUES:
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
                if identifier.lower() not in MISSING_CONFIG_VALUES:
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
        if path.lower() in MISSING_CONFIG_VALUES:
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
    args: object,
    errors: list[str],
    warnings: list[str],
) -> dict[str, set[str]]:
    if args.calibration_regions.strip().lower() in MISSING_CONFIG_VALUES:
        warnings.append(
            "No generated calibration mask supplied; only configured exclusions will filter regions"
        )
    else:
        calibration_regions_kind = genomic_resource_kind(args.calibration_regions)
        if calibration_regions_kind not in {"bed", "tsv"}:
            errors.append(
                "Generated calibration regions must be a BED or TSV file; received "
                f"{args.calibration_regions!r}"
            )
        else:
            validate_indexed_table(
                args.calibration_regions,
                "Calibration regions file",
                errors,
                warnings,
                getattr(args, "calibration_regions_index", None),
            )

    unit_description = f"Configured {args.unit} testing-unit BED file"
    validate_indexed_table(
        args.testing_unit_bed,
        unit_description,
        errors,
        warnings,
        getattr(args, "testing_unit_index", None),
    )
    unit_chromosomes = validate_bed_regions(
        args.testing_unit_bed,
        unit_description,
        errors,
        required_columns=4,
        require_unit_id=True,
    )

    exclusion_paths = optional_config_paths(args.variants_exclude)
    par_paths = optional_config_paths(args.par_regions)
    if not exclusion_paths and not par_paths:
        warnings.append(
            "Neither reference.variants_exclude nor reference.par_regions was "
            "provided; exclusion filtering will be skipped"
        )

    exclusion_indexes = list(getattr(args, "variant_exclude_indexes", []))
    par_indexes = list(getattr(args, "par_region_indexes", []))
    configured_exclusions = exclusion_paths + par_paths
    configured_indexes = exclusion_indexes + par_indexes
    exclusion_chromosomes: set[str] = set()
    for position, path in enumerate(configured_exclusions, start=1):
        index_path = (
            configured_indexes[position - 1]
            if position <= len(configured_indexes)
            else None
        )
        description = f"Exclusion file {position}"
        kind = genomic_resource_kind(path)
        if kind is None:
            errors.append(
                f"{description} must use a .bed, .tsv, .txt, or .vcf extension: {path}"
            )
            continue
        if kind == "vcf":
            validate_vcf_resource(path, description, errors, warnings, index_path)
            continue

        validate_indexed_table(path, description, errors, warnings, index_path)
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
            if gzi_path is None or str(gzi_path).strip().lower() in MISSING_CONFIG_VALUES:
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
        if not cadd_indel or str(cadd_indel).strip().lower() in MISSING_CONFIG_VALUES:
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
