"""Finalize a versioned, provisional null-distribution bundle."""
import csv
import gzip
import json
import math
import os
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_gzip(path, header, rows):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with gzip.open(temporary, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    os.replace(temporary, destination)


def read_bed(path):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                fields = line.rstrip().split("\t")
                chrom, start, end = fields[:3]
                unit = fields[3] if len(fields) > 3 else ""
                rows.append((chrom, int(start), int(end), unit))
    return rows


def main():
    units = read_tsv(snakemake.input.units)
    diagnostics = read_tsv(snakemake.input.diagnostics)
    roles = ("het", "hom")
    by_unit = defaultdict(dict)
    for row in diagnostics:
        by_unit[row["testing_unit_id"]][row["role"]] = row
    excluded, eligible = [], []
    for row in units:
        unit = row["testing_unit_id"]
        if row["status"] != "eligible":
            excluded.append((unit, row["status"]))
            continue
        failed = [role for role in roles if by_unit[unit].get(role, {}).get("status") != "ok"]
        if failed:
            excluded.append((unit, "zero_event_mass_" + "_".join(failed)))
        else:
            eligible.append(unit)
    eligible = sorted(eligible)
    lookups = {}
    for role, path in (("het", snakemake.input.lookup_het), ("hom", snakemake.input.lookup_hom)):
        rows = read_tsv(path)
        grouped = defaultdict(list)
        for row in rows:
            if row["testing_unit_id"] in eligible:
                value = float(row["survival"])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"Invalid {role} survival value")
                grouped[row["testing_unit_id"]].append((int(row["bin"]), value))
        table = []
        for unit in eligible:
            values = [value for _, value in sorted(grouped[unit])]
            if len(values) != int(snakemake.params.bin_count) or any(a < b for a, b in zip(values, values[1:])):
                raise ValueError(f"Invalid {role} lookup grid for {unit}")
            table.append((unit, *values))
        lookups[role] = table
    write_gzip(snakemake.output.lookup_het, ["testing_unit_id", *range(1, int(snakemake.params.bin_count) + 1)], lookups["het"])
    write_gzip(snakemake.output.lookup_hom, ["testing_unit_id", *range(1, int(snakemake.params.bin_count) + 1)], lookups["hom"])
    write_gzip(snakemake.output.excluded_units, ["testing_unit_id", "exclusion_reason"], excluded)
    eligible_bed = [row for row in read_bed(snakemake.input.testing_units) if row[3] in set(eligible)]
    write_gzip(snakemake.output.eligible_bed, ["#CHROM", "START", "END", "testing_unit_id"], eligible_bed)
    par_by_chrom = defaultdict(list)
    for chrom, start, end, _ in read_bed(snakemake.input.par_regions):
        par_by_chrom[chrom.removeprefix("chr").upper()].append((start, end))
    calibration = []
    for chrom, start, end, _ in read_bed(snakemake.input.coverage):
        pieces = [(start, end)]
        for par_start, par_end in par_by_chrom[chrom.removeprefix("chr").upper()]:
            next_pieces = []
            for left, right in pieces:
                if par_end <= left or par_start >= right:
                    next_pieces.append((left, right))
                else:
                    if left < par_start: next_pieces.append((left, par_start))
                    if par_end < right: next_pieces.append((par_end, right))
            pieces = next_pieces
        calibration.extend((chrom, left, right, "coverage_minus_par") for left, right in pieces)
    write_gzip(snakemake.output.calibration_regions, ["#CHROM", "START", "END", "source"], calibration)
    metadata = {"contract": "generated_bundle_v1", "bin_count": int(snakemake.params.bin_count), "bin_width": float(snakemake.params.bin_width), "score_max": float(snakemake.params.score_max), "unobserved_policy": snakemake.params.policy, "mutation_conversion": snakemake.params.conversion if snakemake.params.policy == "mutation_rate" else None, "hemizygote_role": "omitted_par_excluded_initial_scope", "eligible_units": eligible, "grid_status": "configured"}
    Path(snakemake.output.metadata).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    report = {"status": "ok", "eligible_unit_count": len(eligible), "excluded_unit_count": len(excluded), "roles": list(roles)}
    Path(snakemake.output.validation).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    Path(snakemake.output.marker).write_text("validated\n")


if __name__ == "__main__":
    main()
