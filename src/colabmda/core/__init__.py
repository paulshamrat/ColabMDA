"""Core domain models, biological state taxonomy, and project hierarchy."""

from colabmda.core.ions import (
    ION_REGISTRY,
    IonDefinition,
    get_ion,
    parse_ion_spec,
)
from colabmda.core.project import (
    ProjectPaths,
    resolve_project_root,
)
from colabmda.core.taxonomy import (
    SUBSTRATE_PRESETS,
    SystemSpec,
)

__all__ = [
    "IonDefinition",
    "ION_REGISTRY",
    "get_ion",
    "parse_ion_spec",
    "SystemSpec",
    "SUBSTRATE_PRESETS",
    "ProjectPaths",
    "resolve_project_root",
]
