from __future__ import annotations

import os
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

from colabmda.core.taxonomy import SystemSpec


def get_param_path(filename: str) -> Path:
    """Retrieve bundled Amber parameter file from colabmda.amber.params."""
    return Path(str(resources.files("colabmda.amber.params").joinpath(filename)))


def _estimate_ion_pairs(pdb_path: Path, padding_angstroms: float, salt_conc_molar: float) -> int:
    """Estimate the number of monovalent ion pairs to reach target molarity."""
    coords = []
    if not pdb_path.exists():
        return 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    continue
    if not coords:
        return 0
    import numpy as np

    pts = np.array(coords)
    span = pts.max(axis=0) - pts.min(axis=0) + 2.0 * padding_angstroms
    volume_a3 = float(np.prod(span)) * 0.80
    num_pairs = int(round(salt_conc_molar * 0.0006022 * volume_a3))
    return max(0, num_pairs)


def generate_tleap_script(
    spec: SystemSpec,
    input_pdb: Path,
    output_prmtop: Path,
    output_inpcrd: Path,
    output_pdb: Path | None = None,
    water_model: str = "opc",
    protein_ff: str = "ff19SB",
    padding_angstroms: float = 14.0,
    salt_conc_molar: float = 0.15,
    salt_type: str = "NaCl",
) -> str:
    """Generate a clean, self-contained tleap input script for Amber system building.

    Supports:
        - ff19SB / ff14SB protein force fields
        - phosaa19SB for phosphorylated residues (pY171 / PTR)
        - OPC / TIP3P water models with Li-Merz 12-6/12-6-4 ions
        - Hu et al. 2024 quantum parameters for ATP, ADP, GDP, GTP
    """
    lines: list[str] = [
        f"# ColabMDA automated tleap script for {spec.state_tag} ({spec.variant})",
        "",
    ]

    # 1. Atom types for Hu 2024 quantum nucleotide parameters
    is_nucleotide = spec.ligand in {"ATP", "ADP", "GDP", "GTP"}
    if is_nucleotide:
        lines.extend(
            [
                'addAtomTypes {{"O3" "O" "sp2"}}',
                'addAtomTypes {{"O2" "O" "sp2"}}',
                'addAtomTypes {{"O"  "O" "sp2"}}',
                'addAtomTypes {{"OW" "O" "sp3"}}',
                'addAtomTypes {{"OY" "O" "sp3"}}',
                "",
            ]
        )

    # 2. Protein force field
    lines.append(f"source leaprc.protein.{protein_ff}")

    # 3. Phosphorylation parameters (if pY171 / PTR is present)
    if spec.is_phosphorylated:
        if "14" in protein_ff:
            lines.append("source leaprc.phosaa14SB")
        else:
            lines.append("source leaprc.phosaa19SB")

    # 4. Water model and ions
    water_lower = water_model.lower()
    lines.extend(
        [
            f"source leaprc.water.{water_lower}",
            f"loadamberparams frcmod.{water_lower}",
            f"loadamberparams frcmod.ionslm_126_{water_lower}",
            f"loadamberparams frcmod.ionslm_1264_{water_lower}",
            "",
        ]
    )

    # 5. Load nucleotide parameters if applicable
    if is_nucleotide and spec.ligand:
        lig_upper = str(spec.ligand).upper()
        prepi_name = f"{lig_upper}-B3.prepi"
        frcmod_name = "ATP-B3.frcmod"  # Shared B3LYP parameters

        prepi_path = get_param_path(prepi_name)
        frcmod_path = get_param_path(frcmod_name)

        if prepi_path.exists():
            lines.append(f'loadAmberPrep "{prepi_path.resolve()}"')
        if frcmod_path.exists():
            lines.append(f'loadAmberParams "{frcmod_path.resolve()}"')
        lines.append("")

    # 6. Load assembled structure
    abs_pdb = Path(input_pdb).resolve()
    lines.extend(
        [
            f'mol = loadpdb "{abs_pdb}"',
            "charge mol",
            "",
        ]
    )

    # 7. Solvation box
    box_type = "OPCBOX" if water_lower == "opc" else "TIP3PBOX"
    lines.extend(
        [
            f"solvateBox mol {box_type} {padding_angstroms:.1f}",
            "",
        ]
    )

    # 8. Neutralization & Salt
    cation = "Na+" if "na" in salt_type.lower() else "K+"
    anion = "Cl-"
    lines.extend(
        [
            "# Neutralize system charge",
            f"addIonsRand mol {cation} 0",
            f"addIonsRand mol {anion} 0",
            "",
        ]
    )

    if salt_conc_molar > 0:
        num_pairs = _estimate_ion_pairs(abs_pdb, padding_angstroms, salt_conc_molar)
        if num_pairs > 0:
            lines.extend(
                [
                    f"# Add physiological ionic strength ({salt_conc_molar} M {salt_type})",
                    f"addIonsRand mol {cation} {num_pairs} {anion} {num_pairs} 6.0",
                    "",
                ]
            )

    # 9. Save Amber parameters
    abs_top = Path(output_prmtop).resolve()
    abs_crd = Path(output_inpcrd).resolve()
    lines.extend(
        [
            f'saveamberparm mol "{abs_top}" "{abs_crd}"',
        ]
    )
    if output_pdb:
        abs_pdb_out = Path(output_pdb).resolve()
        lines.append(f'savepdb mol "{abs_pdb_out}"')
    lines.extend(
        [
            "charge mol",
            "quit",
            "",
        ]
    )

    return "\n".join(lines)


def run_tleap(
    script_content: str,
    workdir: Path,
    tleap_binary: str | None = None,
    keep_script: bool = False,
) -> tuple[int, str]:
    """Execute tleap in a specific working directory."""
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    exe = tleap_binary or shutil.which("tleap")
    if not exe:
        py_tleap = Path(sys.executable).parent / "tleap"
        if py_tleap.exists():
            exe = str(py_tleap)
    if not exe:
        for env_key in ["AMBERHOME", "AMBER_HOME"]:
            val = os.environ.get(env_key)
            if val and (Path(val) / "bin" / "tleap").exists():
                exe = str(Path(val) / "bin" / "tleap")
                break
    if not exe:
        raise RuntimeError("tleap executable not found in PATH or environment.")

    script_file = workdir / "tleap.in"
    script_file.write_text(script_content, encoding="utf-8")

    cmd = [exe, "-f", "tleap.in"]
    proc = subprocess.run(
        cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    log_text = proc.stdout + "\n" + proc.stderr

    if not keep_script and script_file.exists():
        script_file.unlink()

    return proc.returncode, log_text


def write_prep_provenance(
    spec: SystemSpec,
    dest_json: Path,
    water_model: str = "opc",
    protein_ff: str = "ff19SB",
    padding_angstroms: float = 14.0,
    salt_conc_molar: float = 0.15,
    salt_type: str = "NaCl",
    extra: dict | None = None,
) -> Path:
    """Write machine-readable provenance metadata for system preparation."""
    import json
    from datetime import datetime, timezone

    dest_json = Path(dest_json).resolve()
    dest_json.parent.mkdir(parents=True, exist_ok=True)

    meta = {
        "spec": {
            "protein": spec.protein,
            "state": spec.state_tag,
            "variant": spec.variant,
            "ligand": spec.ligand,
            "ions": spec.ions,
            "phos_residues": spec.phos_residues,
            "peptide_name": spec.peptide_name,
        },
        "parameters": {
            "protein_ff": protein_ff,
            "water_model": water_model,
            "padding_angstroms": padding_angstroms,
            "salt_concentration_molar": salt_conc_molar,
            "salt_type": salt_type,
        },
        "provenance": {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "engine": "ColabMDA",
            "generator": "colabmda.amber.leap",
        },
    }
    if extra:
        meta.update(extra)

    dest_json.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return dest_json
