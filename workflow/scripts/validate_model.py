"""Validate model resources through injected Snakemake inputs and parameters."""

from __future__ import annotations

from types import SimpleNamespace

from apply.validate_inputs import (
    parse_model_settings,
    validate_lookup_files,
    validate_lookup_unit_overlap,
    write_report,
)


def main() -> None:
    args = SimpleNamespace(
        score_max=str(snakemake.params.score_max),
        assembly=str(snakemake.params.assembly),
        unit=str(snakemake.params.unit),
        hemizygote_model=str(snakemake.params.hemizygote_model),
        include_indels=str(snakemake.params.include_indels),
    )
    errors: list[str] = []
    warnings: list[str] = []
    settings = parse_model_settings(args, errors)
    configured = {
        "het": str(snakemake.input.lookup_het),
        "hom": str(snakemake.input.lookup_hom),
        "hem_chrX": str(snakemake.input.lookup_hem_chrx[0]) if snakemake.input.lookup_hem_chrx else "",
    }
    checked = validate_lookup_files(configured, settings, errors)
    overlap = validate_lookup_unit_overlap(snakemake.input.testing_unit_bed, configured, settings, errors)
    write_report({"status": "failed" if errors else "ok", "errors": errors, "warnings": warnings, "summary": {"unit": args.unit, "assembly": args.assembly, "hemizygote_model": settings["hemizygote"], "include_indels": settings["include_indels"], "lookup_files_checked": checked, "lookup_testing_unit_overlap": overlap}}, snakemake.output.report, create_parents=False)


if __name__ == "__main__":
    main()
