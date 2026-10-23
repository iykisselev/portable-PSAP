"""Validate the certified partitioned bundle before cohort work."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from _library import bind
from apply_bindings import (
    ModelValidationInputs,
    ModelValidationParameters,
    ReportOutput,
)
from bundle_certificate import BundleCertificateError, verify_bundle_certificate
from partitioned_bundle_contract import PARTITIONED_NULL_BUNDLE_CONTRACT
from validation_reporting import build_validation_report, write_validation_report


def main() -> None:
    inputs = bind(ModelValidationInputs, snakemake.input)
    outputs = bind(ReportOutput, snakemake.output)
    params = bind(ModelValidationParameters, snakemake.params)
    errors = []
    try:
        component_keys = tuple(params.lookup_component_keys)
        if len(component_keys) != len(inputs.lookup_tables):
            raise ValueError("lookup component keys do not match paths")
        components = dict(zip(component_keys, inputs.lookup_tables, strict=True))
        if len(inputs.resolved_partitions) != len(params.partition_ids) or len(
            inputs.eligible_partitions
        ) != len(params.partition_ids):
            raise ValueError(
                "partition BED components do not match certified partition IDs"
            )
        components.update(
            {
                f"resolved.{partition_id}": path
                for partition_id, path in zip(
                    params.partition_ids, inputs.resolved_partitions, strict=True
                )
            }
        )
        components.update(
            {
                f"eligible.{partition_id}": path
                for partition_id, path in zip(
                    params.partition_ids, inputs.eligible_partitions, strict=True
                )
            }
        )
        components.update(
            {
                "calibration_regions": inputs.lookup_calibration_mask,
                "excluded_alleles": inputs.lookup_excluded_alleles[0],
                "unit_partition_status": inputs.partition_statuses,
                "testing_unit_id_overlap": inputs.testing_unit_id_overlap,
            }
        )
        metadata = verify_bundle_certificate(
            contract=PARTITIONED_NULL_BUNDLE_CONTRACT,
            metadata_path=inputs.lookup_metadata,
            validation_path=inputs.lookup_validation,
            marker_path=inputs.lookup_marker,
            components=components,
        )
        if metadata["assembly"] != params.assembly:
            errors.append("Generated lookup assembly differs from reference.assembly")
        if metadata["testing_unit_type"] != params.unit:
            errors.append("Generated lookup testing-unit type differs from model.unit")
        if bool(params.hemizygote_model) != (
            "hem_chrX" in metadata["model_assumptions"]["roles"]
        ):
            errors.append("Generated lookup HEM model differs from model.hemizygote")
    except (BundleCertificateError, ValueError, KeyError, IndexError) as exc:
        errors.append(f"Could not validate partitioned lookup bundle: {exc}")
    write_validation_report(
        build_validation_report(
            errors=errors,
            warnings=[],
            summary={
                "lookup_contract": params.lookup_contract,
                "lookup_files_checked": len(inputs.lookup_tables),
            },
            sections=("model",),
        ),
        outputs.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
