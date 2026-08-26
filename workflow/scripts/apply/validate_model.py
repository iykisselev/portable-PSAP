#!/usr/bin/env python3
"""Validate model settings and their required lookup resources."""

from __future__ import annotations

import argparse
from pathlib import Path

from validate_inputs import (
    parse_model_settings,
    validate_lookup_files,
    validate_lookup_unit_overlap,
    write_report,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--score-max", required=True)
    result.add_argument("--assembly", required=True)
    result.add_argument("--unit", required=True)
    result.add_argument("--hemizygote-model", required=True)
    result.add_argument("--include-indels", required=True)
    result.add_argument("--lookup-het", required=True)
    result.add_argument("--lookup-hom", required=True)
    result.add_argument("--lookup-hem-chrx", default="")
    result.add_argument("--testing-unit-bed", required=True)
    result.add_argument("--out-json", required=True)
    return result


def main() -> int:
    args = parser().parse_args()
    errors: list[str] = []
    warnings: list[str] = []
    settings = parse_model_settings(args, errors)
    configured_lookups = {
        "het": args.lookup_het,
        "hom": args.lookup_hom,
        "hem_chrX": args.lookup_hem_chrx,
    }
    lookup_paths = validate_lookup_files(
        configured_lookups,
        settings,
        errors,
    )
    overlap = validate_lookup_unit_overlap(
        args.testing_unit_bed,
        configured_lookups,
        settings,
        errors,
    )
    report = {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "unit": args.unit,
            "assembly": args.assembly,
            "hemizygote_model": settings["hemizygote"],
            "include_indels": settings["include_indels"],
            "lookup_files_checked": lookup_paths,
            "lookup_testing_unit_overlap": overlap,
        },
    }
    write_report(report, args.out_json)
    print(
        f"Model validation wrote {Path(args.out_json)} with "
        f"{len(errors)} error(s) and {len(warnings)} warning(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
