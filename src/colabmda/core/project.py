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
        ├── 01_structure/ (or legacy structure/)
        │   ├── protein/
        │   ├── template/
        │   ├── ligand/
        │   ├── peptide/
        │   └── variants/
        ├── 02_simulation/ (or legacy simulation/)
        │   └── <state>/<variant>/
        │       ├── 01_prep/ (or prep/)
        │       ├── 02_eq/   (or eq/)
        │       └── 03_prod/ (or r1/, r2/...)
        └── 03_analysis/ (or legacy analysis/)
            └── <state>/
                ├── <variant>/ (r1/, r2/, aggregate/)
                └── compare/
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
        """Resolve 01_structure or structure."""
        p01 = self.root / "01_structure"
        p_legacy = self.root / "structure"
        if p01.exists() and not p_legacy.exists():
            return p01
        if p_legacy.exists():
            return p_legacy
        return p01 if p01.exists() else p_legacy

    @property
    def structure_01(self) -> Path:
        return self.root / "01_structure"

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
    def variants_dir(self) -> Path:
        return self.structure_dir / "variants"

    @property
    def simulation_base(self) -> Path:
        """Resolve 02_simulation or simulation."""
        p02 = self.root / "02_simulation"
        p_legacy = self.root / "simulation"
        if p02.exists() and not p_legacy.exists():
            return p02
        if p_legacy.exists():
            return p_legacy
        return p02 if p02.exists() else p_legacy

    @property
    def simulation_02(self) -> Path:
        return self.root / "02_simulation"

    @property
    def analysis_base(self) -> Path:
        """Resolve 03_analysis or analysis."""
        p03 = self.root / "03_analysis"
        p_legacy = self.root / "analysis"
        if p03.exists() and not p_legacy.exists():
            return p03
        if p_legacy.exists():
            return p_legacy
        return p03 if p03.exists() else p_legacy

    @property
    def analysis_03(self) -> Path:
        return self.root / "03_analysis"

    def variant_dir(self, spec: SystemSpec) -> Path:
        """Root directory for a specific state and variant."""
        return self.simulation_base / spec.state_tag / spec.variant

    def prep_dir(self, spec: SystemSpec) -> Path:
        """Preparation folder containing model coordinates and Amber topologies."""
        v_dir = self.variant_dir(spec)
        if (v_dir / "01_prep").exists():
            return v_dir / "01_prep"
        if (v_dir / "prep").exists():
            return v_dir / "prep"
        if "02_simulation" in str(self.simulation_base):
            return v_dir / "01_prep"
        return v_dir / "prep"

    def model_pdb(self, spec: SystemSpec) -> Path:
        # Check prep folder first
        p_model = self.prep_dir(spec) / "model.pdb"
        if p_model.exists():
            return p_model
        # Check variants folder: variants/<variant>/model.pdb or variants/<variant>.pdb
        var_dir_model = self.variants_dir / spec.variant / "model.pdb"
        if var_dir_model.exists():
            return var_dir_model
        var_flat = self.variants_dir / f"{spec.variant}.pdb"
        if var_flat.exists():
            return var_flat
        return p_model

    def prepared_prmtop(self, spec: SystemSpec) -> Path:
        return self.prep_dir(spec) / "system.prmtop"

    def prepared_inpcrd(self, spec: SystemSpec) -> Path:
        return self.prep_dir(spec) / "system.inpcrd"

    def prepared_solvated(self, spec: SystemSpec) -> Path:
        """Solvated visualization coordinate PDB (solvated.pdb or system_solvated.pdb)."""
        p_single = self.prep_dir(spec) / "solvated.pdb"
        if p_single.exists():
            return p_single
        p_compat = self.prep_dir(spec) / "system_solvated.pdb"
        if p_compat.exists():
            return p_compat
        return p_single

    def prep_json(self, spec: SystemSpec) -> Path:
        """Provenance and audit metadata for system preparation."""
        return self.prep_dir(spec) / "prep.json"

    def eq_dir(self, spec: SystemSpec) -> Path:
        """Equilibration folder (em, nvt, npt, qc)."""
        v_dir = self.variant_dir(spec)
        if (v_dir / "02_eq").exists():
            return v_dir / "02_eq"
        if (v_dir / "eq").exists():
            return v_dir / "eq"
        if "02_simulation" in str(self.simulation_base):
            return v_dir / "02_eq"
        return v_dir / "eq"

    def eq_pdb(self, spec: SystemSpec) -> Path:
        """Final equilibrated coordinate PDB (eq.pdb or equilibrated.pdb)."""
        p_eq = self.eq_dir(spec) / "eq.pdb"
        return p_eq if p_eq.exists() else (self.eq_dir(spec) / "equilibrated.pdb")

    def prod_dir(self, spec: SystemSpec) -> Path:
        """Production base folder for this variant."""
        v_dir = self.variant_dir(spec)
        p03 = v_dir / "03_prod"
        return p03 if p03.exists() else v_dir

    def replica_dir(self, spec: SystemSpec, replica: str = "r1") -> Path:
        """Production replica folder (pure data: chk, dcd)."""
        norm_rep = f"r{replica}" if str(replica).isdigit() else str(replica).lower()
        v_dir = self.variant_dir(spec)
        p03 = v_dir / "03_prod"
        if p03.exists():
            return p03 / norm_rep
        return v_dir / norm_rep

    def analysis_dir(self, spec: SystemSpec) -> Path:
        """Analysis folder for this specific state and variant."""
        canon = self.analysis_base / spec.state_tag / spec.variant
        if canon.exists():
            return canon
        if "02_simulation" in str(self.simulation_base):
            return canon
        legacy = self.variant_dir(spec) / "analysis"
        if legacy.exists():
            return legacy
        if self.analysis_base.exists():
            return canon
        return legacy

    def compare_dir(self, state_tag: str) -> Path:
        """Cross-variant comparison folder for a state."""
        return self.analysis_base / state_tag / "compare"

    def find_template(self, filename: str) -> Path | None:
        """Search for a template or structure file flexibly across structure subfolders."""
        candidates = [
            self.structure_dir / "template" / filename,
            self.structure_dir / "protein" / filename,
            self.structure_dir / "ligand" / filename,
            self.structure_dir / "peptide" / filename,
            self.structure_dir / "variants" / filename,
            self.structure_dir / filename,
        ]
        # Also search one subfolder level deeper (e.g. template/8FP5/8FP5.pdb or variants/WT/model.pdb)
        for sub in ["template", "protein", "variants", "ligand", "peptide"]:
            subpath = self.structure_dir / sub
            if subpath.is_dir():
                for child in subpath.iterdir():
                    if child.is_dir():
                        candidates.append(child / filename)
        for c in candidates:
            if c.exists():
                return c
        return None

    def init_layout(self, numbered: bool = True) -> None:
        """Create standard top-level project folders."""
        for p in [
            self.scripts_dir,
            self.root / "structure" / "protein",
            self.root / "structure" / "template",
            self.root / "structure" / "ligand",
            self.root / "structure" / "peptide",
            self.root / "structure" / "variants",
            self.root / "simulation",
            self.root / "analysis",
        ]:
            p.mkdir(parents=True, exist_ok=True)

        if numbered:
            for p in [
                self.root / "01_structure" / "protein",
                self.root / "01_structure" / "template",
                self.root / "01_structure" / "ligand",
                self.root / "01_structure" / "peptide",
                self.root / "01_structure" / "variants",
                self.root / "02_simulation",
                self.root / "03_analysis",
            ]:
                p.mkdir(parents=True, exist_ok=True)
