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
    """Standardized, decoupled directory hierarchy for ColabMDA simulation projects.

    Layout:
        <project_root>/
        ├── config.yaml
        ├── scripts/
        │   ├── run.sh
        │   └── submit.slurm
        ├── structure/
        │   ├── protein/
        │   ├── template/
        │   ├── ligand/
        │   └── peptide/
        └── simulation/
            └── <state>/<variant>/
                ├── prep/ (model.pdb, system.prmtop, system.inpcrd)
                ├── eq/   (em, nvt, npt, qc)
                ├── r1/   (chunks, checkpoints, logs)
                └── analysis/
    """

    root: Path

    def __post_init__(self):
        self.root = Path(self.root).resolve()

    @property
    def config_yaml(self) -> Path:
        return self.root / "config.yaml"

    @property
    def scripts_dir(self) -> Path:
        return self.root / "scripts"

    @property
    def structure_dir(self) -> Path:
        return self.root / "structure"

    @property
    def protein_dir(self) -> Path:
        return self.structure_dir / "protein"

    @property
    def template_dir(self) -> Path:
        return self.structure_dir / "template"

    @property
    def ligand_dir(self) -> Path:
        return self.structure_dir / "ligand"

    @property
    def peptide_dir(self) -> Path:
        return self.structure_dir / "peptide"

    @property
    def simulation_base(self) -> Path:
        return self.root / "simulation"

    def variant_dir(self, spec: SystemSpec) -> Path:
        """Root directory for a specific state and variant (e.g. simulation/atp_2mg_eb2/WT)."""
        return self.simulation_base / spec.state_tag / spec.variant

    def prep_dir(self, spec: SystemSpec) -> Path:
        """Preparation folder containing model coordinates and Amber topologies."""
        return self.variant_dir(spec) / "prep"

    def model_pdb(self, spec: SystemSpec) -> Path:
        return self.prep_dir(spec) / "model.pdb"

    def prepared_prmtop(self, spec: SystemSpec) -> Path:
        return self.prep_dir(spec) / "system.prmtop"

    def prepared_inpcrd(self, spec: SystemSpec) -> Path:
        return self.prep_dir(spec) / "system.inpcrd"

    def eq_dir(self, spec: SystemSpec) -> Path:
        """Equilibration folder (em, nvt, npt, qc)."""
        return self.variant_dir(spec) / "eq"

    def replica_dir(self, spec: SystemSpec, replica: str = "r1") -> Path:
        """Production replica folder (pure data: chk, dcd)."""
        norm_rep = f"r{replica}" if str(replica).isdigit() else str(replica).lower()
        return self.variant_dir(spec) / norm_rep

    def analysis_dir(self, spec: SystemSpec) -> Path:
        """Analysis folder for this specific state and variant."""
        return self.variant_dir(spec) / "analysis"

    def find_template(self, filename: str) -> Path | None:
        """Search for a template or structure file flexibly across structure subfolders."""
        candidates = [
            self.structure_dir / "template" / filename,
            self.structure_dir / "protein" / filename,
            self.structure_dir / "ligand" / filename,
            self.structure_dir / "peptide" / filename,
            self.structure_dir / filename,
        ]
        for c in candidates:
            if c.exists():
                return c
        return None

    def init_layout(self) -> None:
        """Create standard top-level project folders."""
        for p in [
            self.scripts_dir,
            self.protein_dir,
            self.template_dir,
            self.ligand_dir,
            self.peptide_dir,
            self.simulation_base,
        ]:
            p.mkdir(parents=True, exist_ok=True)
