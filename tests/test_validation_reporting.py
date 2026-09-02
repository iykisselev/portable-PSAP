"""Contracts for the shared, schema-stable validation report barrier."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from validation_reporting import (  # noqa: E402
    ValidationReportError,
    build_validation_report,
    merge_validation_reports,
    read_validation_report,
    require_validation_success,
    write_validation_report,
)


class ValidationReportingTests(unittest.TestCase):
    def test_reports_merge_with_unique_diagnostics_and_atomic_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            cohort = root / "cohort.json"
            resources = root / "resources.json"
            merged = root / "validation.json"
            marker = root / "validation.ok"
            write_validation_report(
                build_validation_report(
                    errors=[],
                    warnings=["shared warning"],
                    summary={"cohort": 1},
                    sections=("cohort",),
                ),
                cohort,
            )
            write_validation_report(
                build_validation_report(
                    errors=[],
                    warnings=["shared warning", "resource warning"],
                    summary={"resources": 1},
                    sections=("resources",),
                ),
                resources,
            )

            write_validation_report(merge_validation_reports([cohort, resources]), merged)
            report = read_validation_report(merged)
            require_validation_success(
                merged,
                marker,
                failure_message="validation failed",
                include_errors=True,
                success_marker="validated\n",
            )

            self.assertEqual(["shared warning", "resource warning"], report["warnings"])
            self.assertEqual(["cohort", "resources"], report["validation_sections"])
            self.assertEqual("validated\n", marker.read_text(encoding="utf-8"))

    def test_invalid_shape_cannot_create_a_success_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            malformed = Path(temporary_dir) / "malformed.json"
            malformed.write_text('{"status": "ok"}\n', encoding="utf-8")

            with self.assertRaisesRegex(ValidationReportError, "missing required"):
                read_validation_report(malformed)

    def test_failed_report_preserves_actionable_errors(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            report = root / "validation.json"
            write_validation_report(
                build_validation_report(
                    errors=["Cohort sample is absent from the VCF"],
                    warnings=[],
                    summary={},
                    sections=("cohort",),
                ),
                report,
            )

            with self.assertRaisesRegex(RuntimeError, "Cohort sample is absent"):
                require_validation_success(
                    report,
                    root / "validation.ok",
                    failure_message="Input validation failed",
                    include_errors=True,
                    success_marker="validated\n",
                )


if __name__ == "__main__":
    unittest.main()
