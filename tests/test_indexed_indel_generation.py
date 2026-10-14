"""Indexed, boundary-spanning indel generation integration coverage."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

import generation_population
from generation_bindings import (  # noqa: E402
    GenerationWindowInputs,
    GenerationWindowOutputs,
    GenerationWindowParameters,
    WindowStatisticsInputs,
    WindowStatisticsOutputs,
    WindowStatisticsParameters,
)
from generation_window import prepare_generation_window  # noqa: E402
from generation_statistics import calculate_window_statistics  # noqa: E402


@unittest.skipUnless(
    generation_population.pysam, "pysam is required for indexed integration coverage"
)
class IndexedIndelGenerationTests(unittest.TestCase):
    def test_indexed_indel_score_and_rate_cross_a_window_boundary_once(self) -> None:
        pysam = generation_population.pysam
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            window_id = "1_000000000000_000000010000"
            fasta = directory / "reference.fa"
            sequence = list("A" * 20_000)
            sequence[9_999:10_001] = list("CA")
            fasta.write_text(">1\n" + "".join(sequence) + "\n", encoding="utf-8")
            pysam.faidx(str(fasta))
            windows = directory / "windows.json"
            windows.write_text(
                json.dumps(
                    {
                        "schema_version": 2,
                        "windows": [
                            {
                                "window_id": window_id,
                                "chromosome": "1",
                                "start": 0,
                                "end": 10_000,
                                "ordinal": 0,
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            buckets = directory / "buckets.json"
            buckets.write_text(
                json.dumps(
                    {
                        "schema_version": 3,
                        "buckets": [
                            {
                                "bucket_id": "00000",
                                "ordinal": 0,
                                "unit_count": 1,
                                "unit_ids": ["BOUNDARY"],
                                "window_ids": [window_id],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            units = directory / "units.parquet"
            intervals = directory / "intervals.parquet"
            mask = directory / "mask.parquet"
            pl.DataFrame(
                {"testing_unit_id": ["BOUNDARY"], "status": ["eligible"]}
            ).write_parquet(units)
            pl.DataFrame(
                {
                    "testing_unit_id": ["BOUNDARY"],
                    "CHROM": ["1"],
                    "START": [9_999],
                    "END": [10_001],
                }
            ).write_parquet(intervals)
            pl.DataFrame({"CHROM": ["1"], "START": [0], "END": [20_000]}).write_parquet(
                mask
            )

            def indexed_table(name: str, rows: str) -> tuple[Path, Path]:
                plain = directory / f"{name}.tsv"
                compressed = directory / f"{name}.tsv.gz"
                plain.write_text(rows, encoding="utf-8")
                pysam.tabix_compress(str(plain), str(compressed), force=True)
                pysam.tabix_index(
                    str(compressed), seq_col=0, start_col=1, end_col=1, force=True
                )
                return compressed, Path(str(compressed) + ".tbi")

            snv_score, snv_score_index = indexed_table(
                "snv_score",
                "".join(f"chr1\t10000\tC\t{alternate}\t0\t15\n" for alternate in "AGT"),
            )
            indel_score, indel_score_index = indexed_table(
                "indel_score", "chr1\t10000\tCA\tC\t0\t25\n"
            )
            snv_rate, snv_rate_index = indexed_table(
                "snv_rate",
                "".join(
                    f"chr1\t10000\tC\t{alternate}\t0.0001\n" for alternate in "AGT"
                ),
            )
            indel_rate, indel_rate_index = indexed_table(
                "indel_rate", "chr1\t10000\tCA\tC\t0.002\n"
            )
            population_plain = directory / "population.vcf"
            population_gz = directory / "population.vcf.gz"
            population_plain.write_text(
                "##fileformat=VCFv4.2\n"
                '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t10000\t.\tC\tA,G,T\t.\tPASS\tAF=0.1,0.1,0.1\n",
                encoding="utf-8",
            )
            pysam.tabix_compress(str(population_plain), str(population_gz), force=True)
            pysam.tabix_index(str(population_gz), preset="vcf", force=True)

            outputs = GenerationWindowOutputs(
                allele_states=directory / "states.parquet",
                allele_membership=directory / "membership.parquet",
                allele_scores=directory / "scores.parquet",
                candidate_audit=directory / "audit.parquet",
            )
            prepare_generation_window(
                GenerationWindowInputs(
                    windows=windows,
                    units=units,
                    unit_intervals=intervals,
                    effective_mask=mask,
                    population_vcf=population_gz,
                    population_vcf_index=str(population_gz) + ".tbi",
                    score_snv=snv_score,
                    score_snv_index=snv_score_index,
                    score_indel=indel_score,
                    score_indel_index=indel_score_index,
                    fasta=fasta,
                    fasta_fai=str(fasta) + ".fai",
                ),
                outputs,
                GenerationWindowParameters(
                    af_tags={"overall": "AF"},
                    x_af_policy="overall_proxy",
                    hemizygote=False,
                    compression="zstd",
                    include_indels=True,
                    indel_max_reference_length=2,
                    indel_max_alternate_length=2,
                ),
                window_id,
            )
            membership = pl.read_parquet(outputs.allele_membership)
            self.assertEqual(
                membership.filter(pl.col("variant_class") == "indel")
                .select("testing_unit_id", "POS")
                .rows(),
                [("BOUNDARY", 10_000)],
            )

            statistics = directory / "statistics.parquet"
            calculate_window_statistics(
                WindowStatisticsInputs(
                    windows=windows,
                    buckets=buckets,
                    allele_states=outputs.allele_states,
                    allele_membership=outputs.allele_membership,
                    allele_scores=outputs.allele_scores,
                    mutation_rate=snv_rate,
                    mutation_rate_index=snv_rate_index,
                    mutation_rate_indel=indel_rate,
                    mutation_rate_indel_index=indel_rate_index,
                ),
                WindowStatisticsOutputs(statistics=statistics),
                WindowStatisticsParameters(
                    x_af_policy="overall_proxy",
                    unobserved_policy="mutation_rate",
                    x_mutation_rate_policy="sex_averaged_proxy",
                    score_max=70,
                    bin_width=5,
                    roles=("het", "hom"),
                    compression="zstd",
                ),
                window_id,
                "00000",
            )
            result = pl.read_parquet(statistics)
            self.assertGreater(result["snv_event_count"].sum(), 0)
            self.assertGreater(result["indel_event_count"].sum(), 0)


if __name__ == "__main__":
    unittest.main()
