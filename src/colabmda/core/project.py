from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from colabmda.core.taxonomy import SystemSpec

DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/research/projects"
ENV_ROOT = "COLABMDA_ROOT"
ENV_PROJECT = "COLABMDA_PROJECT_DIR"


def resolve_project_root(
    root: str | Path | None = None,
    project_name: str | None = None,
    use_drive: bool = True,
) -> Path:
    """Resolve the base project directory using standardized precedence.

    Priority:
        1. Explicit root path if provided.
        2. Subfolder under root if project_name provided.
        3. Environment variables (COLABMDA_PROJECT_DIR or COLABMDA_ROOT).
        4. Colab Drive mounted path (/content/drive/MyDrive/research/projects/<project>).
        5. Current working directory ($CWD).
    """
    # 1. Explicit root provided
    if root:
        r = Path(root).resolve()
        return r / project_name if project_name else r

    # 2. Check environment variables
    for env_var in [ENV_PROJECT, ENV_ROOT]:
        val = os.environ.get(env_var)
        if val:
            r = Path(val).resolve()
            return r / project_name if project_name else r

    cwd = Path.cwd().resolve()

    # 3. If in Colab temporary /content folder, check Drive
    if str(cwd) == "/content" and use_drive and os.path.exists("/content/drive/MyDrive"):
        base_drive = Path(DEFAULT_DRIVE_ROOT)
        if project_name:
            return base_drive / project_name
        if base_drive.exists() and (base_drive / "research" / "projects").exists():
            return base_drive / "research" / "projects"
        return Path("/content/drive/MyDrive")

    # 4. Standard Current Working Directory
    return cwd / project_name if project_name else cwd


@dataclass
class ProjectPaths:
    """Standardized directory hierarchy for ColabMDA simulation projects."""

    root: Path

    def __post_init__(self):
        self.root = Path(self.root).resolve()

    @property
    def data_dir(self) -> Path:
        return self.root / "data"

    @property
    def raw_dir(self) -> Path:
        return self.data_dir / "raw"

    def models_dir(self, spec: SystemSpec) -> Path:
        return self.data_dir / "models" / spec.state_tag

    def model_pdb(self, spec: SystemSpec) -> Path:
        return self.models_dir(spec) / f"{spec.variant}.pdb"

    def prepared_dir(self, spec: SystemSpec) -> Path:
        return self.data_dir / "prepared" / spec.state_tag / spec.variant

    def prepared_prmtop(self, spec: SystemSpec) -> Path:
        return self.prepared_dir(spec) / "system.prmtop"

    def prepared_inpcrd(self, spec: SystemSpec) -> Path:
        return self.prepared_dir(spec) / "system.inpcrd"

    def simulation_dir(self, spec: SystemSpec, stage_or_replica: str = "equil") -> Path:
        return self.data_dir / "simulation" / spec.state_tag / spec.variant / stage_or_replica

    def processed_dir(self, spec: SystemSpec, stream: str = "keep_ligand") -> Path:
        return self.data_dir / "processed" / spec.state_tag / spec.variant / stream

    def analysis_dir(self, spec: SystemSpec, replica: str | None = None) -> Path:
        base = self.data_dir / "analysis" / spec.state_tag / spec.variant
        return base / replica if replica else base

    def init_layout(self) -> None:
        """Create standard top-level project folders."""
        for p in [self.raw_dir, self.data_dir / "models", self.data_dir / "prepared"]:
            p.mkdir(parents=True, exist_ok=True)
