"""Deployment bindings shared by isolated apply and generation branches."""

import yaml

DEPLOYMENT_DIR = Path(workflow.basedir)
CONDA_ENV_DIR = DEPLOYMENT_DIR / "envs"


def _load_deployment_map(filename):
    with (DEPLOYMENT_DIR / filename).open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


ENVMODULES = _load_deployment_map("envmodules.yaml")
CONTAINERS = _load_deployment_map("containers.yaml")


def mod(*tools):
    """Expand logical tool names into an ordered, de-duplicated module list."""
    modules = []
    for tool in tools:
        for module in ENVMODULES[tool]:
            if module not in modules:
                modules.append(module)
    return modules


def sif(tool):
    """Return the logical tool's pinned Apptainer/Singularity image reference."""
    return CONTAINERS[tool]
