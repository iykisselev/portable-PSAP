"""Resolved lookup and sample-input bindings for apply rules.

Reference paths are supplied by the ``ResourceSpec`` registry in
``resources.smk`` so validation and staging cannot resolve different sources.
"""


INPUT_BINDINGS = APPLY_PLAN.inputs
LOOKUP_SOURCE = APPLY_PLAN.lookup_source
