"""Compatibility facade for the pre-Stage-7 generation helper module.

New code imports focused modules (``generation_partitions``,
``generation_population``, ``generation_scores``, ``generation_reference``,
and ``generation_probability``). These re-exports preserve the historical
public Python API for external users while avoiding a flag-day migration.
"""

try:
    from generation_manifests import PARTITION_MANIFEST_SCHEMA_VERSION as PARTITION_SCHEMA_VERSION
    from generation_partitions import GenomePartition, build_partitions, partition_from_manifest, read_fasta_index, read_partition_manifest
    from generation_population import indexed_lines, open_text, population_info_definition, population_info_number, population_rows, selected_af_values
    from generation_probability import fixed_width_bin_expression, role_probability_expression, unit_bucket_frame
    from generation_reference import ReferenceReader, membership_rows
    from generation_scores import ALLELE_COLUMNS, mutation_rate_rows, score_rows
    from io_contracts import write_json_atomic, write_parquet_atomic
except ModuleNotFoundError:
    from .generation_manifests import PARTITION_MANIFEST_SCHEMA_VERSION as PARTITION_SCHEMA_VERSION
    from .generation_partitions import GenomePartition, build_partitions, partition_from_manifest, read_fasta_index, read_partition_manifest
    from .generation_population import indexed_lines, open_text, population_info_definition, population_info_number, population_rows, selected_af_values
    from .generation_probability import fixed_width_bin_expression, role_probability_expression, unit_bucket_frame
    from .generation_reference import ReferenceReader, membership_rows
    from .generation_scores import ALLELE_COLUMNS, mutation_rate_rows, score_rows
    from .io_contracts import write_json_atomic, write_parquet_atomic


__all__ = [
    "ALLELE_COLUMNS", "GenomePartition", "PARTITION_SCHEMA_VERSION", "ReferenceReader",
    "build_partitions", "fixed_width_bin_expression", "indexed_lines", "membership_rows",
    "mutation_rate_rows", "open_text", "partition_from_manifest", "population_info_definition",
    "population_info_number", "population_rows", "read_fasta_index", "read_partition_manifest",
    "role_probability_expression", "score_rows", "selected_af_values", "unit_bucket_frame",
    "write_json_atomic", "write_parquet_atomic",
]
