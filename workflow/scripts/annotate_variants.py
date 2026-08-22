"""Annotate a filtered VCF with CADD and testing units via injected Snakemake I/O."""

import shlex
import subprocess


def main() -> None:
    include_indels = bool(snakemake.params.include_indels)
    indel_annotation = (
        "bcftools annotate "
        f"--annotations {shlex.quote(str(snakemake.input.cadd_indel[0]))} "
        "--columns CHROM,POS,REF,ALT,-,+INFO/PSAP_CADD_PHRED "
        "--pair-logic exact -Ou -"
        if include_indels
        else "bcftools view -Ou -"
    )
    info_header = (
        "##INFO=<ID=TESTING_UNIT,Number=.,Type=String,"
        f'Description="PSAP {snakemake.params.testing_unit} unit IDs">'
    )
    query_format = (
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/PSAP_CADD_PHRED"
        "\\t%INFO/TESTING_UNIT\\n"
    )
    command = f"""
{{
    bcftools annotate \\
        --annotations {shlex.quote(str(snakemake.input.cadd_snv))} \\
        --columns CHROM,POS,REF,ALT,-,INFO/PSAP_CADD_PHRED \\
        --header-line '##INFO=<ID=PSAP_CADD_PHRED,Number=1,Type=Float,Description="CADD PHRED score">' \\
        --pair-logic exact \\
        -Ou \\
        {shlex.quote(str(snakemake.input.vcf))} | \\
    {indel_annotation} | \\
    bcftools annotate \\
        --annotations {shlex.quote(str(snakemake.input.unit_bed))} \\
        --columns CHROM,FROM,TO,INFO/TESTING_UNIT \\
        --header-line {shlex.quote(info_header)} \\
        --merge-logic TESTING_UNIT:unique \\
        --threads {snakemake.threads} \\
        -Oz \\
        -o {shlex.quote(str(snakemake.output.vcf))} \\
        -

    bcftools index \\
        --threads {snakemake.threads} \\
        --force \\
        {shlex.quote(str(snakemake.output.vcf))}

    {{
        printf 'CHROM\\tPOS\\tREF\\tALT\\tCADD_PHRED\\tunit\\n'
        bcftools query \\
            --format {shlex.quote(query_format)} \\
            {shlex.quote(str(snakemake.output.vcf))}
    }} | bgzip --threads {snakemake.threads} --stdout > {shlex.quote(str(snakemake.output.annotations))}
}}"""
    with open(snakemake.log[0], "wt", encoding="utf-8") as log:
        subprocess.run(command, shell=True, executable="/bin/bash", stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    main()
