"""Class-aware indel event assembly and aggregation invariants."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_event_assembly import (  # noqa: E402
    CLASS_DIAGNOSTIC_SCHEMA,
    EVENT_SCHEMA,
    EventAssemblyError,
    assemble_indel_events,
    event_class_diagnostics,
    role_event_frames,
)
from generation_probability import merge_statistics, statistics_for_role  # noqa: E402


class GenerationEventAssemblyTests(unittest.TestCase):
    @staticmethod
    def inputs() -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
        memberships = pl.DataFrame(
            [
                ("unit-split", "indel", "1", 10, "A", "AT"),
                ("unit-split", "indel", "X", 20, "C", "CT"),
                ("unit-split", "indel", "1", 30, "T", "TA"),
                ("excluded", "indel", "1", 40, "G", "GA"),
            ],
            schema={
                "testing_unit_id": pl.String,
                "variant_class": pl.String,
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
            },
            orient="row",
        )
        scores = pl.DataFrame(
            [
                ("1", 10, "A", "AT", 10.0),
                ("X", 20, "C", "CT", 20.0),
                ("1", 30, "T", "TA", 30.0),
                ("1", 40, "G", "GA", 40.0),
            ],
            schema={"CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String, "score": pl.Float64},
            orient="row",
        )
        states = pl.DataFrame(
            [
                ("1", 10, "A", "AT", "observed", 0.2, None, None),
                ("X", 20, "C", "CT", "observed", 0.3, 0.25, 0.4),
                ("1", 30, "T", "TA", "absent", None, None, None),
                ("1", 40, "G", "GA", "rejected", None, None, None),
            ],
            schema={
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
                "allele_state": pl.String,
                "af_overall": pl.Float64,
                "af_female": pl.Float64,
                "af_male": pl.Float64,
            },
            orient="row",
        )
        assignments = pl.DataFrame(
            [
                ("1", 10, "A", "AT", None, "observed_or_rejected"),
                ("X", 20, "C", "CT", None, "observed_or_rejected"),
                ("1", 30, "T", "TA", 0.01, "absent_mutation_rate"),
                ("1", 40, "G", "GA", None, "observed_or_rejected"),
            ],
            schema={
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
                "mu": pl.Float64,
                "assignment_source": pl.String,
            },
            orient="row",
        )
        return memberships, scores, states, assignments

    def test_exact_indel_assembly_uses_sex_specific_and_fallback_frequencies(self) -> None:
        events = assemble_indel_events(
            *self.inputs(), x_af_policy="sex_specific", x_mutation_rate_policy="sex_averaged_proxy"
        )
        self.assertEqual(events.schema, EVENT_SCHEMA)
        self.assertEqual(events.height, 3)
        by_position = {row["POS"]: row for row in events.iter_rows(named=True)}
        self.assertEqual((by_position[20]["effective_af_standard"], by_position[20]["effective_af_hem"]), (0.25, 0.4))
        self.assertEqual((by_position[30]["effective_af_standard"], by_position[30]["assignment_source"]), (0.01, "absent_mutation_rate"))
        diagnostics = event_class_diagnostics(events, ("het", "hom", "hem_chrX"))
        self.assertEqual(diagnostics.schema, CLASS_DIAGNOSTIC_SCHEMA)
        self.assertEqual(
            {(row["role"], row["variant_class"], row["event_count"]) for row in diagnostics.iter_rows(named=True)},
            {("het", "indel", 3), ("hom", "indel", 3), ("hem_chrX", "indel", 1)},
        )

    def test_absent_x_indel_requires_explicit_mutation_rate_proxy_policy(self) -> None:
        memberships, scores, states, assignments = self.inputs()
        states = pl.concat(
            [states, pl.DataFrame([("X", 21, "A", "AT", "absent", None, None, None)], schema=states.schema, orient="row")]
        )
        memberships = pl.concat(
            [memberships, pl.DataFrame([("unit-split", "indel", "X", 21, "A", "AT")], schema=memberships.schema, orient="row")]
        )
        scores = pl.concat(
            [scores, pl.DataFrame([("X", 21, "A", "AT", 25.0)], schema=scores.schema, orient="row")]
        )
        assignments = pl.concat(
            [assignments, pl.DataFrame([("X", 21, "A", "AT", 0.001, "absent_mutation_rate")], schema=assignments.schema, orient="row")]
        )
        with self.assertRaisesRegex(EventAssemblyError, "configured policy is reject"):
            assemble_indel_events(
                memberships, scores, states, assignments,
                x_af_policy="overall_proxy", x_mutation_rate_policy="reject",
            )

    def test_split_unit_statistics_reduce_to_the_same_one_unit_result(self) -> None:
        events = assemble_indel_events(
            *self.inputs(), x_af_policy="sex_specific", x_mutation_rate_policy="sex_averaged_proxy"
        )
        het = role_event_frames(events, ("het",))["het"].with_columns(pl.lit(0).alias("unit_bucket"))
        direct = statistics_for_role(het, "het", 70, 5).drop("unit_bucket")
        left = statistics_for_role(het.filter(pl.col("POS") <= 20), "het", 70, 5).drop("unit_bucket")
        right = statistics_for_role(het.filter(pl.col("POS") > 20), "het", 70, 5).drop("unit_bucket")
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            direct_path, left_path, right_path = (directory / name for name in ("direct.parquet", "left.parquet", "right.parquet"))
            direct.write_parquet(direct_path)
            left.write_parquet(left_path)
            right.write_parquet(right_path)
            self.assertEqual(merge_statistics([direct_path]), merge_statistics([left_path, right_path]))


if __name__ == "__main__":
    unittest.main()
