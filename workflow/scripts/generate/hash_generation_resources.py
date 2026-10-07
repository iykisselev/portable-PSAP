"""Snakemake adapter for deferred generation-resource provenance hashing."""

import sys

sys.path.insert(0, snakemake.scriptdir)
from _bootstrap import bootstrap

bootstrap(snakemake)

from _library import bind, optional_path

from generation_bindings import GenerationResourceIdentityOutputs
from generation_resource_identities import hash_generation_resources


def main() -> None:
    """Hash every raw input represented in generated-bundle metadata."""
    resources = {
        "population_vcf": str(snakemake.input.population_vcf),
        "score_snv": str(snakemake.input.score_snv),
        "testing_units": str(snakemake.input.testing_units),
        "par_regions": str(snakemake.input.par_regions),
        "fasta": str(snakemake.input.fasta),
    }
    optional_resources = {
        "exclusion": optional_path(snakemake.input.exclusion),
        "mutation_rate": optional_path(snakemake.input.mutation_rate),
        "mutation_rate_indel": optional_path(snakemake.input.mutation_rate_indel),
        "score_indel": optional_path(snakemake.input.score_indel),
    }
    resources.update({name: path for name, path in optional_resources.items() if path})
    resources.update({
        str(identifier): str(path)
        for identifier, path in snakemake.params.exclusions.items()
    })
    hash_generation_resources(
        resources,
        bind(GenerationResourceIdentityOutputs, snakemake.output),
    )


if __name__ == "__main__":
    main()
