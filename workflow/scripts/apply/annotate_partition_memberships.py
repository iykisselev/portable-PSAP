"""Join annotated alleles to certified resolved partition ownership intervals."""

import csv
import gzip
from collections import defaultdict
from pathlib import Path


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def bed_rows(path, partition_id):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError(f"Invalid certified BED row in {path}")
            yield row[0], int(row[1]), int(row[2]), row[3], partition_id


def main():
    partition_ids = tuple(snakemake.params.partition_ids)
    if len(partition_ids) != len(snakemake.input.resolved) or len(partition_ids) != len(snakemake.input.eligible):
        raise ValueError("Certified partition BED inputs do not match configured order")
    ownership = defaultdict(list)
    eligible = set()
    statuses = {}
    for order, (partition_id, path) in enumerate(zip(partition_ids, snakemake.input.resolved, strict=True)):
        for chrom, start, end, unit, _ in bed_rows(path, partition_id):
            ownership[(chrom, unit)].append((start, end, order, partition_id))
    for partition_id, path in zip(partition_ids, snakemake.input.eligible, strict=True):
        for chrom, _start, _end, unit, _ in bed_rows(path, partition_id):
            eligible.add((chrom, unit, partition_id))
    for row in rows(snakemake.input.statuses):
        key = (row["CHROM"], row["testing_unit_id"], row["partition_id"])
        if key in statuses:
            raise ValueError(f"Duplicate certified partition status for {key}")
        statuses[key] = (row["status"], row.get("failed_roles_json", "[]"))
    for values in ownership.values():
        values.sort(key=lambda value: (value[0], value[1], value[2]))
    destination = Path(str(snakemake.output.memberships))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(destination, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "CHROM", "POS", "REF", "ALT", "CADD_PHRED", "unit", "partition_id",
            "partition_status", "failed_roles_json",
        ])
        for allele in rows(snakemake.input.annotations):
            chrom, position, ref = allele["CHROM"], int(allele["POS"]), allele["REF"]
            start, end = position - 1, position - 1 + len(ref)
            for (candidate_chrom, unit), intervals in ownership.items():
                if candidate_chrom != chrom:
                    continue
                matches = [item for item in intervals if item[0] < end and item[1] > start]
                if not matches:
                    continue
                _, _, _order, partition_id = min(matches, key=lambda item: item[2])
                status_key = (chrom, unit, partition_id)
                if status_key not in statuses:
                    raise ValueError(f"Missing certified partition status for {status_key}")
                status, failed_roles = statuses[status_key]
                # Keep ownership evidence even for an ineligible pair.  Apply
                # scoring filters it later, while cohort reporting uses it to
                # retain an explicit unavailable-partition row in G.
                writer.writerow([
                    chrom, position, ref, allele["ALT"], allele["CADD_PHRED"], unit,
                    partition_id, status, failed_roles,
                ])


if __name__ == "__main__":
    main()
