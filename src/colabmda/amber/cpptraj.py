from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def generate_cpptraj_strip_script(
    topology: Path,
    trajectories: list[Path],
    output_traj: Path,
    output_top: Path | None = None,
    mode: str = "keep_ligand",
    autoimage: bool = True,
) -> str:
    """Generate a cpptraj script for autoimaging and trajectory stripping.

    Modes:
        - 'keep_ligand': Keeps protein, peptide, nucleotide, and active-site ions.
                         Strips bulk water and salt (:WAT,Na+,Cl-,K+).
        - 'ca_backbone': Strips everything except protein C-alpha atoms (@CA)
                         for cross-state comparative PCA and RMSD with apo.
        - 'protein_only': Strips water, salt, and all ligands/ions.
    """
    lines: list[str] = [
        f'parm "{Path(topology).resolve()}"',
    ]

    for traj in trajectories:
        lines.append(f'trajin "{Path(traj).resolve()}"')

    if autoimage:
        lines.append("autoimage")

    if mode == "keep_ligand":
        # Strip solvent and neutralizing ions, keeping all solutes
        lines.append("strip :WAT,Na+,Cl-,K+,CIO,IB+ outprefix stripped")
    elif mode == "ca_backbone":
        # Keep only C-alpha atoms
        lines.append("strip !@CA outprefix ca_only")
    elif mode == "protein_only":
        lines.append("strip !(:1-9999 & !:WAT,Na+,Cl-,K+,ATP,ADP,GDP,MG,MN)")

    abs_out_traj = Path(output_traj).resolve()
    lines.append(f'trajout "{abs_out_traj}"')

    if output_top:
        abs_out_top = Path(output_top).resolve()
        lines.append(f'parmwrite out "{abs_out_top}"')

    lines.append("run")
    lines.append("quit")
    return "\n".join(lines)


def run_cpptraj(
    script_content: str,
    workdir: Path,
    cpptraj_binary: str | None = None,
) -> tuple[int, str]:
    """Execute cpptraj with the given script content."""
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    exe = cpptraj_binary or shutil.which("cpptraj")
    if not exe:
        py_cpptraj = Path(sys.executable).parent / "cpptraj"
        if py_cpptraj.exists():
            exe = str(py_cpptraj)
    if not exe:
        for env_key in ["AMBERHOME", "AMBER_HOME"]:
            val = os.environ.get(env_key)
            if val and (Path(val) / "bin" / "cpptraj").exists():
                exe = str(Path(val) / "bin" / "cpptraj")
                break
    if not exe:
        raise RuntimeError("cpptraj executable not found in PATH or environment.")

    with tempfile.TemporaryDirectory(prefix="colabmda-cpptraj-") as sandbox:
        script_file = Path(sandbox) / "cpptraj.in"
        script_file.write_text(script_content, encoding="utf-8")
        proc = subprocess.run(
            [exe, "-i", str(script_file)], cwd=sandbox, capture_output=True, text=True
        )
    return proc.returncode, proc.stdout + "\n" + proc.stderr
