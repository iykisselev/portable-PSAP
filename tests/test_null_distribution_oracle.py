"""Independent fixtures for the planned null-distribution generator.

This oracle intentionally does not import production workflow code. It
enumerates Bernoulli outcomes for tiny units so a later probability recurrence
cannot validate itself by sharing helpers or expected values.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from itertools import product
import json
from pathlib import Path
import unittest


FIXTURE_DIR = Path(__file__).parent / "null_distribution"
getcontext().prec = 50


def load_json(name: str) -> dict:
    with (FIXTURE_DIR / name).open(encoding="utf-8") as handle:
        return json.load(handle)


def enumerate_maximum(events: list[dict[str, str | int]]) -> tuple[Decimal, dict[int, Decimal]]:
    """Enumerate every event subset and return no-event and maximum-score mass."""
    no_event = Decimal("0")
    maximum_mass: dict[int, Decimal] = {}
    probabilities = [Decimal(str(event["q"])) for event in events]
    scores = [int(event["score"]) for event in events]

    for outcome in product((False, True), repeat=len(events)):
        probability = Decimal("1")
        observed_scores: list[int] = []
        for happened, q, score in zip(outcome, probabilities, scores):
            probability *= q if happened else Decimal("1") - q
            if happened:
                observed_scores.append(score)
        if not observed_scores:
            no_event += probability
        else:
            maximum = max(observed_scores)
            maximum_mass[maximum] = maximum_mass.get(maximum, Decimal("0")) + probability
    return no_event, maximum_mass


def parse_info(text: str) -> dict[str, str]:
    if text == ".":
        return {}
    return dict(item.split("=", 1) for item in text.split(";"))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def overlap_length(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    total = 0
    for left_start, left_end in merge_intervals(left):
        for right_start, right_end in merge_intervals(right):
            total += max(0, min(left_end, right_end) - max(left_start, right_start))
    return total


class NullDistributionOracleTests(unittest.TestCase):
    def test_exhaustive_event_enumeration_matches_hand_computed_cases(self) -> None:
        fixture = load_json("probability_oracle.json")
        for case in fixture["cases"]:
            with self.subTest(case=case["name"]):
                no_event, maximum_mass = enumerate_maximum(case["events"])
                expected = case["expected"]
                self.assertEqual(no_event, Decimal(expected["none"]))
                self.assertEqual(
                    maximum_mass,
                    {int(score): Decimal(mass) for score, mass in expected["max_mass"].items()},
                )
                event_mass = Decimal("1") - no_event
                for threshold, survival in expected["conditional_survival"].items():
                    tail = sum(mass for score, mass in maximum_mass.items() if score >= int(threshold))
                    self.assertLessEqual(
                        abs(tail / event_mass - Decimal(survival)),
                        Decimal("1e-27"),
                    )

    def test_legacy_and_provisional_generated_bins_are_deliberately_distinct(self) -> None:
        fixture = load_json("probability_oracle.json")
        score_max = Decimal(str(fixture["score_max"]))
        legacy = fixture["bin_mapping_examples"]["external_legacy_v0"]
        generated = fixture["bin_mapping_examples"]["generated_bundle_v1"]

        legacy_grid = [
            score_max * Decimal(index) / Decimal(legacy["column_count"] - 1)
            for index in range(legacy["column_count"])
        ]
        generated_width = score_max / Decimal(generated["bin_count"])

        for score_text, expected_column in legacy["expected_columns"].items():
            score = min(max(Decimal(score_text), Decimal("0")), score_max)
            legacy_column = max(index for index, edge in enumerate(legacy_grid, start=1) if edge <= score)
            self.assertEqual(legacy_column, expected_column)
        for score_text, expected_column in generated["expected_columns"].items():
            score = min(max(Decimal(score_text), Decimal("0")), score_max)
            generated_column = min(int(score // generated_width) + 1, generated["bin_count"])
            self.assertEqual(generated_column, expected_column)

        self.assertEqual(legacy["expected_columns"]["5"], 1)
        self.assertEqual(generated["expected_columns"]["5"], 2)

    def test_population_panel_fixture_preserves_alt_frequency_identity(self) -> None:
        expected = load_json("fixture_expectations.json")["af_assignment"]
        records = []
        with (FIXTURE_DIR / "population_af.vcf").open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    records.append(line.rstrip("\n").split("\t"))

        def extract(tag: str) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
            observed: dict[str, str] = {}
            invalid: dict[str, str] = {}
            excluded: dict[str, str] = {}
            for chrom, position, _, ref, alts, _, filter_value, info_text in records:
                info = parse_info(info_text)
                alternatives = alts.split(",")
                if filter_value != "PASS":
                    for alt in alternatives:
                        excluded[f"{chrom}:{position}:{ref}:{alt}"] = "non_pass_filter"
                    continue
                value = info.get(tag)
                for index, alt in enumerate(alternatives):
                    key = f"{chrom}:{position}:{ref}:{alt}"
                    if value is None or value == ".":
                        invalid[key] = "missing_configured_af"
                    elif len(alternatives) > 1 and "," not in value:
                        invalid[key] = "scalar_af_on_multiallelic_record"
                    else:
                        observed[key] = value.split(",")[index]
            return observed, invalid, excluded

        observed, invalid, excluded = extract(expected["configured_tag"])

        self.assertEqual(observed, expected["observed"])
        self.assertEqual(invalid, expected["invalid"])
        self.assertEqual(excluded, expected["excluded"])
        self.assertNotIn(expected["absent_eligible_key"], observed)
        _, scalar_invalid, _ = extract("AF_SCALAR")
        self.assertEqual(
            {key: reason for key, reason in scalar_invalid.items() if key.startswith("chr1:50:")},
            expected["scalar_tag_invalid"],
        )

    def test_interval_fixture_uses_nonpar_denominator_and_retained_bases(self) -> None:
        expected = load_json("fixture_expectations.json")["retention"]
        units: dict[str, list[tuple[str, int, int]]] = {}
        exclusions: dict[str, list[tuple[int, int]]] = {}
        with (FIXTURE_DIR / "testing_units.bed").open(encoding="utf-8") as handle:
            for line in handle:
                chrom, start, end, unit = line.rstrip("\n").split("\t")
                units.setdefault(unit, []).append((chrom, int(start), int(end)))
        with (FIXTURE_DIR / "exclusion.bed").open(encoding="utf-8") as handle:
            for line in handle:
                chrom, start, end = line.rstrip("\n").split("\t")
                exclusions.setdefault(chrom, []).append((int(start), int(end)))
        with (FIXTURE_DIR / "par_regions.bed").open(encoding="utf-8") as handle:
            par = [tuple(line.rstrip("\n").split("\t")) for line in handle]

        for unit, intervals in units.items():
            by_chrom: dict[str, list[tuple[int, int]]] = {}
            for chrom, start, end in intervals:
                by_chrom.setdefault(chrom, []).append((start, end))
            self.assertEqual(len(by_chrom), 1)
            chrom, unit_intervals = next(iter(by_chrom.items()))
            unit_bases = {position for start, end in unit_intervals for position in range(start, end)}
            par_bases = {position for par_chrom, start, end in par if par_chrom == chrom for position in range(int(start), int(end))}
            excluded_bases = {position for start, end in exclusions.get(chrom, []) for position in range(start, end)}
            supported = unit_bases - par_bases
            retained = supported - excluded_bases
            self.assertEqual(len(supported), expected[unit]["supported_bases"])
            self.assertEqual(len(retained), expected[unit]["retained_bases"])
            if "fraction" in expected[unit]:
                self.assertEqual(Decimal(len(retained)) / Decimal(len(supported)), Decimal(expected[unit]["fraction"]))
            if "par_union_bases" in expected[unit]:
                par_intervals = [(int(start), int(end)) for par_chrom, start, end in par if par_chrom == chrom]
                self.assertEqual(overlap_length(unit_intervals, par_intervals), expected[unit]["par_union_bases"])


if __name__ == "__main__":
    unittest.main()
