"""Contracts for bounded normalized indel catalog candidates."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
import warnings

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_indels import (  # noqa: E402
    INDEL_MEMBERSHIP_SCHEMA,
    INDEL_MUTATION_ASSIGNMENT_SCHEMA,
    INDEL_MUTATION_RATE_SCHEMA,
    INDEL_REJECTION_SCHEMA,
    INDEL_SCORE_SCHEMA,
    IndelCandidateError,
    indel_candidate_memberships,
    indel_mutation_assignments,
    indel_mutation_rate_rows,
    indel_population_states,
    indel_score_rows,
)
from generation_population import population_rows  # noqa: E402


class InMemoryReference:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        if chromosome != "1" or start < 0 or end > len(self.sequence):
            raise ValueError("outside test reference")
        return self.sequence[start:end]


class GenerationIndelTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.resource = self.directory / "indels.tsv"
        self.reference = InMemoryReference("CAGTAC")

    def rows(self, content: str, *, start: int = 0, end: int = 6) -> pl.DataFrame:
        self.resource.write_text(content, encoding="utf-8")
        return indel_score_rows(
            self.resource,
            None,
            "1",
            start,
            end,
            self.reference,
            max_reference_length=3,
            max_alternate_length=3,
        )

    def test_reader_uses_anchor_owner_window_and_preserves_duplicates(self) -> None:
        rows = self.rows(
            "chr1\t2\tA\tAT\t0\t15\n"
            "1\t2\tA\tAT\t0\t15\n"
            "1\t7\tC\tCA\t0\t21\n"
        )
        self.assertEqual(rows.schema, INDEL_SCORE_SCHEMA)
        self.assertEqual(rows.height, 2)
        self.assertEqual(rows["POS"].to_list(), [2, 2])
        self.assertEqual(rows["SPAN_START"].to_list(), [1, 1])
        self.assertEqual(rows["SPAN_END"].to_list(), [2, 2])

    def test_reader_rejects_invalid_catalog_rows(self) -> None:
        cases = {
            "not_indel": ("1\t2\tA\tG\t0\t15\n", "not a length-changing"),
            "symbolic": ("1\t2\tA\t<DEL>\t0\t15\n", "SNV or length-changing"),
            "not_minimal": ("1\t2\tAGT\tAT\t0\t15\n", "not minimal"),
            "not_left_aligned": ("1\t2\tA\tAA\t0\t15\n", "not left-aligned"),
            "reference_mismatch": ("1\t2\tC\tCT\t0\t15\n", "does not match"),
            "length": ("1\t2\tAGTA\tA\t0\t15\n", "exceeds configured"),
            "score": ("1\t2\tA\tAT\t0\tnan\n", "non-finite"),
        }
        for name, (content, message) in cases.items():
            with self.subTest(name=name), self.assertRaisesRegex(IndelCandidateError, message):
                self.rows(content)

    def test_full_span_mask_and_unit_membership_are_distinct(self) -> None:
        candidates = pl.DataFrame(
            [
                ("1", 10, "AAA", "A", 9, 12, 20.0),
                ("1", 7, "AA", "A", 6, 8, 21.0),
                ("1", 9, "A", "AT", 8, 9, 22.0),
            ],
            schema=INDEL_SCORE_SCHEMA,
            orient="row",
        )
        memberships, rejected = indel_candidate_memberships(
            candidates,
            mask=[(8, 13)],
            unit_intervals=[("unit-a", 10, 11), ("unit-b", 11, 13)],
        )
        self.assertEqual(memberships.schema, INDEL_MEMBERSHIP_SCHEMA)
        self.assertEqual(
            [(row["testing_unit_id"], row["POS"], row["variant_class"])
             for row in memberships.sort("testing_unit_id").iter_rows(named=True)],
            [("unit-a", 10, "indel"), ("unit-b", 10, "indel")],
        )
        self.assertEqual(rejected.schema, INDEL_REJECTION_SCHEMA)
        self.assertEqual(
            {(row["POS"], row["rejection_reason"]) for row in rejected.iter_rows(named=True)},
            {
                (7, "span_not_contained_in_calibration_mask"),
                (9, "no_eligible_testing_unit_overlap"),
            },
        )

    def test_population_state_join_keeps_zero_imputation_and_absence_distinct(self) -> None:
        candidates = pl.DataFrame(
            [
                ("1", 2, "A", "AT", 1, 2, 20.0),
                ("1", 3, "G", "GA", 2, 3, 21.0),
                ("1", 4, "T", "TA", 3, 4, 22.0),
            ],
            schema=INDEL_SCORE_SCHEMA,
            orient="row",
        )
        population = self.directory / "population.vcf"
        population.write_text(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=AF_NFE,Number=A,Type=Float,Description=\"AF\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\t.\tA\tAT\t.\tPASS\tAF_NFE=0.2\n"
            "1\t3\t.\tG\tGA\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            panel = population_rows(population, None, "1", 0, 6, {"overall": "AF_NFE"})
        states = indel_population_states(candidates, panel)
        by_position = {row["POS"]: row for row in states.iter_rows(named=True)}
        self.assertEqual((by_position[2]["allele_state"], by_position[2]["af"]), ("observed", 0.2))
        self.assertEqual(
            (by_position[3]["state_reason"], by_position[3]["af_zero_imputed_tags"]),
            ("selected_af_zero_imputed", "AF_NFE"),
        )
        self.assertEqual(
            (by_position[4]["allele_state"], by_position[4]["state_reason"]),
            ("absent", "not_present_in_population_panel"),
        )
        self.assertTrue(any("using AF=0" in str(warning.message) for warning in caught))

    def test_indel_mutation_rate_fallback_is_exact_and_class_specific(self) -> None:
        states = pl.DataFrame(
            [
                ("1", 2, "A", "AT", "observed"),
                ("1", 4, "T", "TA", "absent"),
                ("1", 5, "A", "AT", "rejected"),
            ],
            schema={"CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String, "allele_state": pl.String},
            orient="row",
        )
        rates_path = self.directory / "indel_rates.tsv"
        rates_path.write_text("chr1\t4\tT\tTA\t0.0001\n", encoding="utf-8")
        rates = indel_mutation_rate_rows(
            rates_path,
            None,
            "1",
            0,
            6,
            self.reference,
            max_reference_length=3,
            max_alternate_length=3,
        )
        self.assertEqual(rates.schema, INDEL_MUTATION_RATE_SCHEMA)
        assigned = indel_mutation_assignments(states, rates, "mutation_rate")
        self.assertEqual(assigned.schema, INDEL_MUTATION_ASSIGNMENT_SCHEMA)
        by_position = {row["POS"]: row for row in assigned.iter_rows(named=True)}
        self.assertEqual((by_position[4]["mu"], by_position[4]["assignment_source"]), (0.0001, "absent_mutation_rate"))
        self.assertIsNone(by_position[2]["mu"])
        zero = indel_mutation_assignments(states, None, "zero")
        self.assertEqual(zero.filter(pl.col("POS") == 4)["mu"].to_list(), [0.0])
        with self.assertRaisesRegex(IndelCandidateError, "missing 1 absent"):
            indel_mutation_assignments(states, pl.DataFrame(schema=INDEL_MUTATION_RATE_SCHEMA), "mutation_rate")

    def test_mutation_rate_reader_rejects_non_indel_or_invalid_rate(self) -> None:
        for content, message in (
            ("1\t4\tT\tC\t0.1\n", "not an indel"),
            ("1\t4\tT\tTA\tnan\n", "invalid mu"),
        ):
            with self.subTest(content=content), self.assertRaisesRegex(IndelCandidateError, message):
                self.resource.write_text(content, encoding="utf-8")
                indel_mutation_rate_rows(
                    self.resource,
                    None,
                    "1",
                    0,
                    6,
                    self.reference,
                    max_reference_length=3,
                    max_alternate_length=3,
                )


if __name__ == "__main__":
    unittest.main()
