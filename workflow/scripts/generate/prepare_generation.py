"""Prepare panel AF, coverage eligibility, and SNV score-universe diagnostics."""

from pathlib import Path
import csv
import json
import sys

# Snakemake adds only this nested script directory to sys.path.  The explicit
# scripts root keeps nested wrappers independent of their relative depth.
SCRIPTS_ROOT = Path(snakemake.scriptdir).parent
sys.path.insert(0, str(SCRIPTS_ROOT))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_prep import (
    allele_key,
    coverage_rows,
    key_text,
    read_bed_units,
    read_fasta,
    read_intervals,
    read_population_af,
    read_scores,
)


def write_tsv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    observed, rejected = read_population_af(snakemake.input.population_vcf, snakemake.params.af_tag)
    units = read_bed_units(snakemake.input.testing_units)
    coverage = read_intervals(snakemake.input.coverage)
    unit_rows = coverage_rows(units, coverage, float(snakemake.params.coverage_minimum))
    sequences = read_fasta(snakemake.input.fasta)
    scores = read_scores(snakemake.input.score_snv)

    observed_rows = [(*key, value, snakemake.params.af_tag) for key, value in sorted(observed.items())]
    rejected_rows = [(*key, reason) for key, reason in sorted(rejected.items())]
    possible, score_rows, score_gaps = set(), [], []
    eligible_units = {unit for unit, _, _, _, _, status in unit_rows if status == "eligible"}
    for unit in sorted(eligible_units):
        for chrom, start, end in units[unit]:
            sequence = sequences.get(chrom)
            if sequence is None:
                raise ValueError(f"Testing-unit chromosome {chrom!r} is absent from generation FASTA")
            if end > len(sequence):
                raise ValueError(f"Testing-unit interval for {unit!r} exceeds generation FASTA")
            for offset in range(start, end):
                ref = sequence[offset]
                if ref not in "ACGT":
                    continue
                for alt in "ACGT":
                    if alt == ref:
                        continue
                    key = allele_key(chrom, offset + 1, ref, alt)
                    if key in possible:
                        continue
                    possible.add(key)
                    score = scores.get(key)
                    if score is None:
                        score_gaps.append((*key, "missing_score"))
                    else:
                        score_rows.append((unit, *key, score))

    write_tsv(snakemake.output.observed_af, ["CHROM", "POS", "REF", "ALT", "af", "af_source_tag"], observed_rows)
    write_tsv(snakemake.output.excluded_alleles, ["CHROM", "POS", "REF", "ALT", "exclusion_reason"], rejected_rows)
    write_tsv(snakemake.output.units, ["testing_unit_id", "CHROM", "union_bases", "covered_bases", "coverage_fraction", "status"], unit_rows)
    write_tsv(snakemake.output.score_universe, ["testing_unit_id", "CHROM", "POS", "REF", "ALT", "score"], score_rows)
    write_tsv(snakemake.output.score_gaps, ["CHROM", "POS", "REF", "ALT", "reason"], score_gaps)
    report = {
        "status": "ok",
        "af_source_tag": snakemake.params.af_tag,
        "observed_af_count": len(observed),
        "rejected_panel_allele_count": len(rejected),
        "eligible_unit_count": len(eligible_units),
        "score_universe_count": len(possible),
        "scored_allele_count": len(score_rows),
        "missing_score_count": len(score_gaps),
        "missing_score_policy": "reported_pending_D4_D6",
    }
    Path(snakemake.output.report).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
