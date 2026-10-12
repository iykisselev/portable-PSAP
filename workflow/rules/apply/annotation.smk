import shlex

ANNOTATED_VCF = ANNOTATION_DIR / f"{TESTING_UNIT}.vcf.gz"
ANNOTATED_VCF_CSI = Path(f"{ANNOTATED_VCF}.csi")
RAW_UNIT_ANNOTATIONS = ANNOTATION_DIR / f"{TESTING_UNIT}.annotated.tsv.gz"
ANNOTATION_QC_REPORT = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.json"
ANNOTATION_QC_OK = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.ok"


rule annotate_variants:
    """
    Attach CADD and the configured testing-unit annotation to retained variants.
    """
    input:
        vcf=FILTERED_VCF,
        vcf_csi=FILTERED_VCF_INDEX,
        cadd_snv=apply_resource("cadd_snv").staged_data,
        cadd_snv_index=apply_resource("cadd_snv").staged_index,
        cadd_indel=APPLY_RESOURCE_REGISTRY.inputs("cadd_indel", "staged_data"),
        cadd_indel_index=APPLY_RESOURCE_REGISTRY.inputs("cadd_indel", "staged_index"),
        unit_bed=apply_resource("testing_units").staged_data,
        unit_bed_index=apply_resource("testing_units").staged_index,
    output:
        vcf=ensure(str(ANNOTATED_VCF), non_empty=True),
        csi=ensure(str(ANNOTATED_VCF_CSI), non_empty=True),
        annotations=ensure(str(RAW_UNIT_ANNOTATIONS), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "annotation" / "annotate_variants.log"),
    benchmark:
        str(BENCHMARK_DIR / "annotate_variants.tsv")
    conda:
        env("bioinformatics")
    container:
        sif("bioinformatics")
    envmodules:
        *mod("bioinformatics"),
    # bcftools annotate processes one stream per process.  Split the indexed,
    # filtered VCF into fixed-width regions and run up to this many independent
    # streams below; bcftools' own --threads option only helps compression.
    threads: 4
    resources:
        mem_mb=policy_resources("annotate_variants").mem_mb,
        runtime=policy_resources("annotate_variants").runtime,
    params:
        partition_size_bp=50_000_000,
        cadd_header=(
            "##INFO=<ID=PSAP_CADD_PHRED,Number=1,Type=Float,"
            'Description="CADD PHRED score">'
        ),
        unit_header=(
            "##INFO=<ID=TESTING_UNIT,Number=.,Type=String,"
            f'Description="PSAP {TESTING_UNIT} unit IDs">'
        ),
        indel_stage=lambda wc, input: (
            "bcftools annotate "
            f"--annotations {shlex.quote(str(input.cadd_indel[0]))} "
            "--columns CHROM,POS,REF,ALT,-,+INFO/PSAP_CADD_PHRED "
            "--pair-logic exact -Ou - |"
            if input.cadd_indel
            else ""
        ),
    shell:
        """
        {{
            set -euo pipefail
            parts_dir=$(mktemp -d "${{TMPDIR:-/tmp}}/psap-annotation.XXXXXX")
            trap 'rm -rf "$parts_dir"' EXIT
            regions="$parts_dir/regions.tsv"

            # Index statistics are in VCF-header contig order.  Emit only
            # regions carrying records, split into fixed-width, non-overlapping
            # one-based intervals so concat restores the original sort order.
            bcftools index --stats {input.vcf:q} |
                awk -v width={params.partition_size_bp} '\
                    $2 > 0 && $3 > 0 {{
                        for (start = 1; start <= $2; start += width) {{
                            end = start + width - 1
                            if (end > $2) end = $2
                            printf "%06d\\t%s:%d-%d\\n", ++ordinal, $1, start, end
                        }}
                    }}' > "$regions"

            annotate_region() {{
                local region="$1"
                local destination="$2"
                if [[ -n "$region" ]]; then
                    # Partition ownership is the VCF anchor POS, not an
                    # indel's reference span.  That makes adjacent regions
                    # disjoint even for boundary-spanning indels.
                    bcftools view --regions "$region" --regions-overlap 0 -Ou {input.vcf:q}
                else
                    # An empty filtered VCF still needs a valid annotated
                    # header and an empty, indexed final VCF.
                    bcftools view -Ou {input.vcf:q}
                fi |
                bcftools annotate \
                    --annotations {input.cadd_snv:q} \
                    --columns CHROM,POS,REF,ALT,-,INFO/PSAP_CADD_PHRED \
                    --header-line {params.cadd_header:q} \
                    --pair-logic exact \
                    -Ou - |
                {params.indel_stage}
                bcftools annotate \
                    --annotations {input.unit_bed:q} \
                    --columns CHROM,FROM,TO,INFO/TESTING_UNIT \
                    --header-line {params.unit_header:q} \
                    --merge-logic TESTING_UNIT:unique \
                    -Ob -o "$destination" -
                bcftools index --force "$destination"
            }}

            declare -a workers=()
            while IFS=$'\\t' read -r ordinal region; do
                annotate_region "$region" "$parts_dir/$ordinal.bcf" &
                workers+=("$!")
                if (( ${{#workers[@]}} == {threads} )); then
                    wait "${{workers[0]}}"
                    workers=("${{workers[@]:1}}")
                fi
            done < "$regions"

            if (( ${{#workers[@]}} == 0 )) && [[ ! -s "$regions" ]]; then
                annotate_region "" "$parts_dir/000001.bcf"
            fi
            for worker in "${{workers[@]}}"; do
                wait "$worker"
            done

            mapfile -t part_files < <(find "$parts_dir" -maxdepth 1 -name '*.bcf' -print | sort)
            bcftools concat --threads {threads} -Oz -o {output.vcf:q} "${{part_files[@]}}"
            bcftools index --threads {threads} --force {output.vcf:q}

            {{
                printf 'CHROM\\tPOS\\tREF\\tALT\\tCADD_PHRED\\tunit\\n'
                bcftools query \
                    --format '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/PSAP_CADD_PHRED\\t%INFO/TESTING_UNIT\\n' \
                    {output.vcf:q}
            }} | bgzip --threads {threads} --stdout > {output.annotations:q}
        }} > {log:q} 2>&1
        """


rule summarize_annotations:
    """Report complete annotation misses before any per-sample scoring."""
    input:
        annotations=rules.annotate_variants.output.annotations,
    output:
        report=ensure(str(ANNOTATION_QC_REPORT), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "annotation" / "summarize_annotations.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("summarize_annotations").mem_mb,
        runtime=policy_resources("summarize_annotations").runtime,
    params:
        include_indels=WORKFLOW_PLAN.model.include_indels,
    script:
        str(APPLY_SCRIPT_DIR / "summarize_annotations.py")


rule require_valid_annotations:
    """Stop the workflow when an enabled annotation source matched nothing."""
    input:
        report=rules.summarize_annotations.output.report,
    output:
        ok=ensure(str(ANNOTATION_QC_OK), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "annotation" / "require_valid_annotations.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("require_valid_annotations").mem_mb,
        runtime=policy_resources("require_valid_annotations").runtime,
    params:
        failure_message="Annotation validation failed",
        include_errors=False,
        success_marker="annotation validation passed\n",
    script:
        str(APPLY_SCRIPT_DIR / "require_valid_report.py")
