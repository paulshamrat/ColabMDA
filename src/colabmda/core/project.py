from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

from colabmda.core.taxonomy import SystemSpec

DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/projects"
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
        4. Colab Drive mounted path (/content/drive/MyDrive/projects/<project>).
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
        if os.path.exists("/content/drive/MyDrive/projects"):
            base = Path("/content/drive/MyDrive/projects")
            return base / project_name if project_name else base
        if os.path.exists("/content/drive/MyDrive/ColabMDA"):
            base = Path("/content/drive/MyDrive/ColabMDA")
            return base / project_name if project_name else base
        if os.path.exists("/content/drive/MyDrive/research/projects"):
            base = Path("/content/drive/MyDrive/research/projects")
            return base / project_name if project_name else base
        base_drive = Path(DEFAULT_DRIVE_ROOT)
        return base_drive / project_name if project_name else base_drive

    # 4. Standard Current Working Directory
    return cwd / project_name if project_name else cwd


@dataclass
class ProjectPaths:
    """Standardized, decoupled directory hierarchy for ColabMDA simulation projects.

    New state-aware CLI runs use data/<state>/structure and
    data/<state>/simulation/<run>/<variant>/{01_prep,02_eq,03_md} and
    data/<state>/analysis/<run>/<variant>/{r1,r2,r3,aggregate}.
    Existing numbered and legacy runs remain readable in their original locations.
    """

    root: Path
    state: str | None = None

    def __post_init__(self):
        self.root = Path(self.root).resolve()
        if self.state and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", self.state):
            raise ValueError("Invalid state")
        self.run_id = os.environ.get("COLABMDA_RUN_ID", "")
        if self.run_id and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", self.run_id):
            raise ValueError("Invalid COLABMDA_RUN_ID")

    @property
    def state_dir(self) -> Path | None:
        """Use state-first paths for new CLI runs; preserve existing run locations."""
        if self.state is None:
            return None
        if self.run_id and any(
            (self.root / "data" / folder / self.run_id).exists()
            for folder in ["02_simulation", "02_pilot"]
        ):
            return None
        if not self.run_id and any(
            (base / self.state).exists()
            for base in [
                self.root / "data/simulation",
                self.root / "simulation",
                self.root / "02_simulation",
            ]
        ):
            return None
        name = "atpmg" if self.state in {"atp_mg", "atpmg", "atp_1mg"} else self.state
        return self.root / "data" / name

    @property
    def config_yaml(self) -> Path:
        return self.root / "config.yaml"

    @property
    def scripts_dir(self) -> Path:
        return self.root / "scripts"

    @property
    def structure_dir(self) -> Path:
        """Resolve data/structure, 01_structure, or structure."""
        if self.state_dir is not None:
            return self.state_dir / "structure"
        numbered = self.root / "data/01_structure"
        if self.run_id or numbered.exists():
            return numbered
        p_data = self.root / "data" / "structure"
        p01 = self.root / "01_structure"
        p_legacy = self.root / "structure"
        if p_data.exists():
            return p_data
        if p01.exists() and not p_legacy.exists():
            return p01
        if p_legacy.exists():
            return p_legacy
        return p_legacy

    @property
    def structure_01(self) -> Path:
        return self.root / "01_structure"

    @property
    def protein_dir(self) -> Path:
        return self.structure_dir / ("models/WT" if self.state_dir is not None else "protein")

    @property
    def template_dir(self) -> Path:
        return self.structure_dir / ("raw" if self.state_dir is not None else "template")

    @property
    def ligand_dir(self) -> Path:
        return self.structure_dir / "ligand"

    @property
    def peptide_dir(self) -> Path:
        return self.structure_dir / "peptide"

    @property
    def variants_dir(self) -> Path:
        return self.structure_dir / ("models" if self.state_dir is not None else "variants")

    @property
    def simulation_base(self) -> Path:
        """Resolve data/simulation, 02_simulation, or simulation."""
        if self.state_dir is not None:
            base = self.state_dir / "simulation"
            return base / self.run_id if self.run_id else base
        if self.run_id:
            numbered = self.root / "data/02_simulation" / self.run_id
            if numbered.exists() or (self.root / "data/00_setup/example.json").exists():
                return numbered
            return self.root / "data/02_pilot" / self.run_id
        p_data = self.root / "data" / "simulation"
        p02 = self.root / "02_simulation"
        p_legacy = self.root / "simulation"
        if p_data.exists():
            return p_data
        if p02.exists() and not p_legacy.exists():
            return p02
        if p_legacy.exists():
            return p_legacy
        return p_legacy

    @property
    def simulation_02(self) -> Path:
        return self.root / "02_simulation"

    @property
    def analysis_base(self) -> Path:
        """Resolve data/analysis, 03_analysis, or analysis."""
        if self.state_dir is not None:
            return self.state_dir / "analysis" / self.run_id
        p_data = self.root / "data" / "analysis"
        p03 = self.root / "03_analysis"
        p_legacy = self.root / "analysis"
        if p_data.exists():
            return p_data
        if p03.exists() and not p_legacy.exists():
            return p03
        if p_legacy.exists():
            return p_legacy
        return p_legacy

    @property
    def analysis_03(self) -> Path:
        return self.root / "03_analysis"

    def variant_dir(self, spec: SystemSpec) -> Path:
        """Root directory for a specific state and variant."""
        if self.state_dir is not None:
            if not self.run_id:
                raise ValueError(
                    "Set COLABMDA_RUN_ID to a fresh run name, or its original name to resume"
                )
            return self.simulation_base / spec.variant
        return self.simulation_base / spec.state_tag / spec.variant

    def prep_dir(self, spec: SystemSpec) -> Path:
        """Preparation folder containing model coordinates and Amber topologies."""
        v_dir = self.variant_dir(spec)
        if self.run_id or self.state_dir is not None:
            return v_dir / "01_prep"
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
        if self.run_id or self.state_dir is not None:
            return v_dir / "02_eq"
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
        if self.run_id or self.state_dir is not None:
            return v_dir / "03_md"
        p03 = v_dir / "03_prod"
        return p03 if p03.exists() else v_dir

    def replica_dir(self, spec: SystemSpec, replica: str = "r1") -> Path:
        norm_rep = f"r{replica}" if str(replica).isdigit() else str(replica).lower()
        if not re.fullmatch(r"r[0-9]+", norm_rep):
            raise ValueError("Replica must be r followed by an integer")
        return self.prod_dir(spec) / norm_rep

    def analysis_dir(self, spec: SystemSpec) -> Path:
        """Analysis folder for this specific state and variant."""
        if self.state_dir is not None:
            self.variant_dir(spec)  # Require the same explicit run ID as simulation.
            return self.analysis_base / spec.variant
        if self.run_id:
            if "02_simulation" in self.simulation_base.parts:
                return self.variant_dir(spec) / "04_analysis"
            return self.variant_dir(spec) / "04_check"
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

    def replica_analysis_dir(self, spec: SystemSpec, replica: str = "r1") -> Path:
        """Canonical per-replica reports; preserve pre-state layouts verbatim."""
        normalized = self.replica_dir(spec, replica).name
        return (
            self.analysis_dir(spec) / normalized
            if self.state_dir is not None
            else self.analysis_dir(spec)
        )

    def migrate_analysis(self, spec: SystemSpec) -> dict | None:
        """Explicitly migrate completed legacy high-level r1 reports, hash verified.

        Never called implicitly by path resolution. Historical numbered examples
        remain unchanged. A receipt preserves the old/new paths and file hashes.
        """
        import hashlib
        import json
        import shutil

        if self.state_dir is None:
            raise ValueError("Analysis migration applies only to state-first runs")
        source = self.variant_dir(spec) / "04_analysis"
        target = self.replica_analysis_dir(spec, "r1")
        if not source.exists():
            return None
        required = [
            "ANALYSIS_DONE",
            "analysis.json",
            "rmsd.csv",
            "rg.csv",
            "rmsf.csv",
            "figures/rmsd_vs_time.png",
            "figures/rg_vs_time.png",
            "figures/rmsf_per_residue.png",
        ]
        if any(
            not (source / name).is_file() or not (source / name).stat().st_size for name in required
        ):
            raise ValueError("Only completed legacy r1 analysis can be migrated")
        if any((source / name).exists() for name in ["r1", "r2", "r3", "aggregate"]):
            raise ValueError(
                "Legacy analysis contains multiple replicas; explicit mapping required"
            )

        def inventory(folder):
            result = {}
            for path in [folder, *folder.rglob("*")]:
                if path.is_symlink():
                    raise ValueError("Refusing symlink in analysis migration")
                if path.is_file():
                    with path.open("rb") as handle:
                        digest = hashlib.sha256()
                        for block in iter(lambda: handle.read(1024 * 1024), b""):
                            digest.update(block)
                    result[str(path.relative_to(folder))] = digest.hexdigest()
            return result

        before = inventory(source)
        if any(Path(n).suffix in {".dcd", ".xtc", ".trr", ".nc", ".chk"} for n in before):
            raise ValueError("Analysis migration refuses heavy simulation data")
        # Refuse parent symlinks and all conflicting files before copying anything.
        for parent in [source, *source.parents, target, *target.parents]:
            if parent.is_symlink():
                raise ValueError("Refusing symlink in analysis destination")
        existing = inventory(target) if target.exists() else {}
        if any(name not in before or digest != before[name] for name, digest in existing.items()):
            raise ValueError("Conflicting analysis destination preserved")
        target.mkdir(parents=True, exist_ok=True)
        markers = [n for n in before if n.endswith("_DONE")]
        for name in markers:
            (target / name).unlink(missing_ok=True)
        for name in before.keys() - set(markers):
            dest = target / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / name, dest)
        if inventory(source) != before or inventory(target) != {
            n: h for n, h in before.items() if n not in markers
        }:
            raise ValueError("Analysis changed during migration; original preserved")
        record = {
            "source": str(source.relative_to(self.root)),
            "destination": str(target.relative_to(self.root)),
            "replica": "r1",
            "files": before,
            "provenance": "Original metadata retained byte-for-byte; paths in historical receipts remain original",
        }
        receipt = (
            self.root
            / "data/00_setup/analysis_migrations"
            / f"{self.state_dir.name}_{self.run_id}_{spec.variant}.json"
        )
        receipt.parent.mkdir(parents=True, exist_ok=True)
        if receipt.exists() and json.loads(receipt.read_text()) != record:
            raise ValueError("Conflicting migration receipt preserved")
        receipt.write_text(json.dumps(record, indent=2) + "\n")
        for name in markers:
            dest = target / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / name, dest)
        if inventory(source) != before or inventory(target) != before:
            for name in markers:
                (target / name).unlink(missing_ok=True)
            raise ValueError("Final migration verification failed; original preserved")
        shutil.rmtree(source)
        return record

    def log_dir(self, spec: SystemSpec) -> Path:
        """Execution logs mirror state/run/variant, independently of scientific data."""
        if self.state_dir is not None:
            return self.root / "logs" / self.state_dir.name / self.run_id / spec.variant
        return self.root / "logs" / "02_simulation" / self.run_id / spec.variant

    def compare_dir(self, state_tag: str) -> Path:
        """Cross-variant comparison folder for a state."""
        if self.state_dir is not None:
            return self.analysis_base / "compare"
        return self.analysis_base / state_tag / "compare"

    def find_template(self, filename: str) -> Path | None:
        """Search for a template or structure file flexibly across structure subfolders."""
        bases = [self.structure_dir]
        if self.state_dir is None:
            bases += [self.root / "data" / "structure", self.root / "structure"]
        candidates = []
        for b in bases:
            if not b.exists():
                continue
            candidates.extend(
                [
                    b / "template" / filename,
                    b / "protein" / filename,
                    b / "ligand" / filename,
                    b / "peptide" / filename,
                    b / "variants" / filename,
                    b / filename,
                ]
            )
            # Also search subfolders (e.g. 4bgq/wt, variants/WT, template/8FP5)
            for sub in ["template", "protein", "variants", "ligand", "peptide", "4bgq"]:
                subpath = b / sub
                if subpath.is_dir():
                    for child in subpath.iterdir():
                        if child.is_dir():
                            candidates.append(child / filename)
            # Comprehensive fallback search across structure directory
            for match in b.rglob(filename):
                if match.is_file():
                    candidates.append(match)
        for c in candidates:
            if c.exists():
                return c
        return None

    def init_layout(self, numbered: bool = True) -> None:
        """Create only the resolved layout; never mirror compatibility trees."""
        if self.state_dir is not None:
            for folder in ["raw", "mapping", "models"]:
                (self.structure_dir / folder).mkdir(parents=True, exist_ok=True)
            self.simulation_base.mkdir(parents=True, exist_ok=True)
            return
        for path in [
            self.scripts_dir,
            self.protein_dir,
            self.template_dir,
            self.ligand_dir,
            self.peptide_dir,
            self.variants_dir,
            self.simulation_base,
        ]:
            path.mkdir(parents=True, exist_ok=True)
        if not self.run_id:
            self.analysis_base.mkdir(parents=True, exist_ok=True)
