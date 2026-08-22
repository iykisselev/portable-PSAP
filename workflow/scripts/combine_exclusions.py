"""Project configured exclusions through injected Snakemake I/O."""

from pathlib import Path

from apply.combine_exclusions import collect_exclusions, write_projected_vcf


def main() -> None:
    intervals, variants = collect_exclusions(snakemake.params.sources)
    write_projected_vcf(Path(snakemake.input.cohort), Path(snakemake.output.vcf), intervals, variants)


if __name__ == "__main__":
    main()
