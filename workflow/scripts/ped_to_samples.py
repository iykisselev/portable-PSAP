#!/usr/bin/env python3
"""Convert the first six PED columns to the Stage 1 sample-sheet contract."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ped", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with Path(args.ped).open("rt", encoding="utf-8", newline="") as source, output.open(
        "wt", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(["family_id", "pid", "sex", "phenotype"])
        for line_number, row in enumerate(csv.reader(source, delimiter="\t"), start=1):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) < 6:
                raise ValueError(
                    f"PED line {line_number} has {len(row)} columns; at least 6 are required"
                )
            writer.writerow([row[0], row[1], row[4], row[5]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
