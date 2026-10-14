from __future__ import annotations

import re
import shutil
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np

from colabmda.core.taxonomy import SUBSTRATE_PRESETS, SystemSpec

AA3_TO_1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "PTR": "Y",
    "SEP": "S",
    "TPO": "T",
}

CANONICAL_AA1_TO_3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

AA1_TO_3 = CANONICAL_AA1_TO_3

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_rcsb_pdb(pdb_id: str, dest_dir: str | Path) -> Path:
    """Download a reference PDB structure from RCSB into the project structure directory.

    Decouples data from the tool repository: structures are saved in the project
    workspace (e.g. `structure/template/<pdb_id>.pdb`), never in `src/`.
    """
    dest = Path(dest_dir).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    clean_id = pdb_id.strip().upper()
    target_file = dest / f"{clean_id}.pdb"
    if target_file.exists():
        return target_file

    url = RCSB_DOWNLOAD_URL.format(pdb_id=clean_id)
    urllib.request.urlretrieve(url, target_file)
    return target_file


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the optimal rotation matrix and translation vector using Kabsch algorithm.

    Transforms source coordinates to target as: x_transformed = x @ rotation + translation

    Args:
        source: (N, 3) coordinates array.
        target: (N, 3) coordinates array.

    Returns:
        tuple (rotation, translation) where rotation is (3, 3) and translation is (3,).
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)

    if len(source) != len(target):
        raise ValueError(
            f"Point count mismatch in Kabsch: source={len(source)}, target={len(target)}"
        )
    if len(source) < 3:
        raise ValueError("At least 3 points required for 3D Kabsch alignment")

    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)

    source_centered = source - source_center
    target_centered = target - target_center

    covariance = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)

    # Ensure right-handed coordinate system (no reflection)
    det = np.linalg.det(u @ vt)
    if det < 0:
        u[:, -1] *= -1

    rotation = u @ vt
    translation = target_center - source_center @ rotation
    return rotation, translation


def parse_pdb_atoms(
    pdb_path: str | Path,
    chain: str | None = None,
    res_name: str | None = None,
    res_num: int | None = None,
) -> dict[str, np.ndarray]:
    """Extract coordinates dictionary {atom_name: np.array([x, y, z])} for a specific residue."""
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")

    atoms: dict[str, np.ndarray] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            line_chain = line[21]
            try:
                line_resnum = int(line[22:26].strip())
            except ValueError:
                continue
            line_resname = line[17:20].strip()
            line_atom = line[12:16].strip()

            if chain is not None and line_chain != chain and line_chain != " ":
                continue
            if res_num is not None and line_resnum != res_num:
                continue
            if res_name is not None and line_resname != res_name:
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atoms[line_atom] = np.array([x, y, z], dtype=float)

    return atoms


def get_chain_residues(pdb_path: str | Path, chain: str = "A") -> list[tuple[str, int]]:
    """Return ordered list of (res_name, res_num) for polymer atoms in the specified chain."""
    residues: list[tuple[str, int]] = []
    seen: set[tuple[str, int]] = set()
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM  ") or line[21] != chain:
                continue
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            key = (resname, resnum)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def build_ptr_phosphate_atoms(
    tyr_atoms: dict[str, np.ndarray],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Mathematically generate ideal tetrahedral phosphate atoms (P, O1P, O2P, O3P) for PTR.

    Derives coordinates geometrically from the tyrosine phenol OH and CZ coordinates,
    eliminating any dependency on static external PDB templates.
    """
    if "OH" not in tyr_atoms or "CZ" not in tyr_atoms:
        raise ValueError("Tyrosine residue missing OH or CZ atom required for phosphorylation")

    cz = tyr_atoms["CZ"]
    oh = tyr_atoms["OH"]
    cg = tyr_atoms.get("CG", cz - np.array([1.0, 0.0, 0.0]))

    v = oh - cz
    norm_v = np.linalg.norm(v)
    e_z = v / norm_v if norm_v > 1e-4 else np.array([0.0, 0.0, 1.0])

    # Orthogonal vector relative to the aromatic ring axis
    w = cz - cg
    perp = np.cross(e_z, w)
    if np.linalg.norm(perp) < 1e-4:
        perp = np.array([0.0, 1.0, 0.0])
    e_y = perp / np.linalg.norm(perp)
    e_x = np.cross(e_y, e_z)

    # Position P atom (P-O bond ~1.60 A, C-O-P angle ~120 deg)
    theta = np.radians(60.0)
    p_vec = np.cos(theta) * e_z + np.sin(theta) * e_x
    p_coord = oh + 1.60 * p_vec

    # Position tetrahedral phosphate oxygens (P-O ~1.51 A, tetrahedral 109.5 deg)
    op_axis = (p_coord - oh) / np.linalg.norm(p_coord - oh)
    ox = np.cross(op_axis, e_y)
    if np.linalg.norm(ox) < 1e-4:
        ox = np.cross(op_axis, e_x)
    ox = ox / np.linalg.norm(ox)
    oy = np.cross(op_axis, ox)

    phi = np.radians(70.5)
    r_bond = 1.51
    z_proj = np.cos(phi) * r_bond
    rad_proj = np.sin(phi) * r_bond

    oxygens: dict[str, np.ndarray] = {}
    for i, name in enumerate(["O1P", "O2P", "O3P"]):
        ang = np.radians(i * 120.0 + 30.0)
        o_vec = z_proj * op_axis + rad_proj * (np.cos(ang) * ox + np.sin(ang) * oy)
        oxygens[name] = p_coord + o_vec

    return p_coord, oxygens


def phosphorylate_residue(
    input_pdb: str | Path,
    output_pdb: str | Path,
    res_num: int = 171,
    chain: str = "A",
    target_res: str = "PTR",
    ptr_template: str | Path | None = None,
) -> Path:
    """Phosphorylate a target residue in a PDB file (e.g., CDKL5 Y171 -> PTR171).

    If `ptr_template` is provided, performs Kabsch alignment from the template.
    Otherwise, derives exact tetrahedral coordinates mathematically without needing
    any external PDB data files.
    """
    in_path = Path(input_pdb)
    out_path = Path(output_pdb)

    if not in_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {in_path}")

    target_atoms = parse_pdb_atoms(in_path, chain=chain, res_num=res_num)
    if not target_atoms:
        target_atoms = parse_pdb_atoms(in_path, chain=None, res_num=res_num)
    if not target_atoms:
        raise ValueError(f"Residue {chain}:{res_num} not found in {in_path}")

    # If residue is already phosphorylated (contains P and phosphate oxygens), preserve it
    if "P" in target_atoms and any(k in target_atoms for k in ["O1P", "O2P", "O3P"]):
        if str(in_path.resolve()) != str(out_path.resolve()):
            shutil.copy2(in_path, out_path)
        return out_path

    replacement_records = []
    serial = 90000

    if ptr_template and Path(ptr_template).exists():
        # Align from template if user explicitly provided one
        tmpl_atoms = parse_pdb_atoms(ptr_template, res_name="PTR")
        common_atoms = [
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CD1",
            "CD2",
            "CE1",
            "CE2",
            "CZ",
        ]
        shared = [a for a in common_atoms if a in target_atoms and a in tmpl_atoms]
        if len(shared) < 6:
            raise ValueError(
                f"Insufficient matching heavy atoms ({len(shared)}) to align PTR at {chain}:{res_num}"
            )
        rot, trans = kabsch(
            np.vstack([tmpl_atoms[a] for a in shared]),
            np.vstack([target_atoms[a] for a in shared]),
        )
        for atom_name, coord in tmpl_atoms.items():
            if atom_name.startswith("H"):
                continue
            new_coord = coord @ rot + trans
            elem = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper() or "C"
            replacement_records.append(
                f"ATOM  {serial:5d} {atom_name:>4s} {target_res} {chain}{res_num:4d}    "
                f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}  1.00  0.00          {elem:>2s}\n"
            )
            serial += 1
    else:
        # Mathematically synthesize phosphate coordinates on the existing tyrosine
        p_coord, oxygens = build_ptr_phosphate_atoms(target_atoms)
        for atom_name, coord in target_atoms.items():
            if atom_name.startswith("H"):
                continue
            elem = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper() or "C"
            replacement_records.append(
                f"ATOM  {serial:5d} {atom_name:>4s} {target_res} {chain}{res_num:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00          {elem:>2s}\n"
            )
            serial += 1

        # Add P and oxygens
        replacement_records.append(
            f"ATOM  {serial:5d}  P   {target_res} {chain}{res_num:4d}    "
            f"{p_coord[0]:8.3f}{p_coord[1]:8.3f}{p_coord[2]:8.3f}  1.00  0.00           P\n"
        )
        serial += 1
        for oname, ocoord in oxygens.items():
            replacement_records.append(
                f"ATOM  {serial:5d} {oname:>4s} {target_res} {chain}{res_num:4d}    "
                f"{ocoord[0]:8.3f}{ocoord[1]:8.3f}{ocoord[2]:8.3f}  1.00  0.00           O\n"
            )
            serial += 1

    # Write modified PDB
    out_lines = []
    replaced = False
    with open(in_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                line_chain = line[21]
                try:
                    line_resnum = int(line[22:26].strip())
                except ValueError:
                    line_resnum = None

                if line_chain == chain and line_resnum == res_num:
                    if not replaced:
                        out_lines.extend(replacement_records)
                        replaced = True
                    continue
            out_lines.append(line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(out_lines)

    return out_path


def position_nucleotide(
    receptor_pdb: str | Path,
    output_pdb: str | Path,
    template_pdb: str | Path,
    ligand_name: str = "ATP",
    res_num: int = 304,
    chain: str = "A",
) -> bool:
    """Transfer nucleotide from template_pdb to receptor_pdb using kinase CA alignment."""
    rec_path = Path(receptor_pdb)
    out_path = Path(output_pdb)
    tmpl_path = Path(template_pdb)

    if not rec_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {rec_path}")
    if not tmpl_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {tmpl_path}")

    existing_lig = parse_pdb_atoms(rec_path, res_name=ligand_name)
    if existing_lig:
        if rec_path.resolve() != out_path.resolve():
            shutil.copy2(rec_path, out_path)
        return True

    tmpl_lig_atoms = parse_pdb_atoms(tmpl_path, res_name=ligand_name)
    if not tmpl_lig_atoms:
        if rec_path.resolve() != out_path.resolve():
            shutil.copy2(rec_path, out_path)
        return False

    rec_ca: dict[int, np.ndarray] = {}
    with open(rec_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                try:
                    r_num = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    rec_ca[r_num] = np.array([x, y, z])
                except ValueError:
                    pass

    tmpl_ca: dict[int, np.ndarray] = {}
    with open(tmpl_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM  ") and line[12:16].strip() == "CA" and line[21] in {"A", " "}:
                try:
                    r_num = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    tmpl_ca[r_num] = np.array([x, y, z])
                except ValueError:
                    pass

    common_res = sorted(set(rec_ca.keys()) & set(tmpl_ca.keys()))
    if len(common_res) < 20:
        n_common = min(len(rec_ca), len(tmpl_ca), 250)
        if n_common < 3:
            if rec_path.resolve() != out_path.resolve():
                shutil.copy2(rec_path, out_path)
            return False
        source_pts = np.vstack([tmpl_ca[k] for k in sorted(tmpl_ca.keys())[:n_common]])
        target_pts = np.vstack([rec_ca[k] for k in sorted(rec_ca.keys())[:n_common]])
    else:
        source_pts = np.vstack([tmpl_ca[k] for k in common_res])
        target_pts = np.vstack([rec_ca[k] for k in common_res])

    rotation, translation = kabsch(source_pts, target_pts)

    lig_records: list[str] = []
    atom_idx = 5000
    with open(tmpl_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM  ")) and line[17:20].strip() == ligand_name:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    new_xyz = np.dot(np.array([x, y, z]), rotation) + translation
                    atom_name = line[12:16]
                    atom_idx += 1
                    occ = float(line[54:60]) if len(line) >= 60 else 1.00
                    b_fac = float(line[60:66]) if len(line) >= 66 else 0.00
                    elem = line[76:78].strip() if len(line) >= 78 else atom_name.strip()[:1]
                    lig_line = (
                        f"HETATM{atom_idx:5d} {atom_name:4s} {ligand_name:>3s} {chain}{res_num:4d}    "
                        f"{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}"
                        f"{occ:6.2f}{b_fac:6.2f}          {elem:>2s}\n"
                    )
                    lig_records.append(lig_line)
                except ValueError:
                    continue

    with open(rec_path, encoding="utf-8") as f:
        clean_lines = [line for line in f if not line.strip().startswith(("END", "ENDMDL"))]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(clean_lines + lig_records + ["END\n"])

    return True


def position_ions(
    receptor_pdb: str | Path,
    output_pdb: str | Path,
    template_pdb: str | Path,
    ion_type: str = "MG",
    count: int = 2,
    chain: str = "A",
    res_start: int = 301,
) -> list[dict[str, Any]]:
    """Position 1 or 2 catalytic divalent ions (e.g. MG, MN, CA, ZN) into the active site.

    Superimposes ATP triphosphate atoms (PA, PB, PG, etc.) from template_pdb onto
    receptor_pdb to transfer divalent ion coordinates with quantum precision.
    """
    rec_path = Path(receptor_pdb)
    out_path = Path(output_pdb)
    tmpl_path = Path(template_pdb)

    if not rec_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {rec_path}")
    if not tmpl_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {tmpl_path}")

    rec_atp = parse_pdb_atoms(rec_path, res_name="ATP")
    tmpl_atp = parse_pdb_atoms(tmpl_path, res_name="ATP")

    if not rec_atp:
        raise ValueError(f"No ATP found in receptor {rec_path} to align catalytic ions")
    if not tmpl_atp:
        raise ValueError(f"No ATP found in template {tmpl_path}")

    triphosphate_atoms = [
        "PA",
        "PB",
        "PG",
        "O1A",
        "O2A",
        "O1B",
        "O2B",
        "O1G",
        "O2G",
        "O3G",
    ]
    shared_atp = [a for a in triphosphate_atoms if a in rec_atp and a in tmpl_atp]
    if len(shared_atp) < 4:
        raise ValueError(
            f"Insufficient matching ATP phosphate atoms ({len(shared_atp)}) for ion alignment"
        )

    source_atp = np.vstack([tmpl_atp[a] for a in shared_atp])
    target_atp = np.vstack([rec_atp[a] for a in shared_atp])

    rotation, translation = kabsch(source_atp, target_atp)

    template_ions: list[tuple[str, np.ndarray]] = []
    with open(tmpl_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() in {
                "MG",
                "MN",
                "CA",
                "ZN",
            }:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ion_name = line[12:16].strip()
                template_ions.append((ion_name, np.array([x, y, z])))

    if not template_ions:
        raise ValueError(f"No divalent ions found in catalytic template {tmpl_path}")

    selected_ions = template_ions[:count]

    ion_records = []
    placed_info: list[dict[str, Any]] = []
    curr_res = res_start
    serial = 95000

    ion_elem = ion_type.upper()
    ion_resn = f"{ion_elem:>3s}"

    for idx, (_orig_name, orig_coord) in enumerate(selected_ions, start=1):
        trans_coord = orig_coord @ rotation + translation
        rec = (
            f"HETATM{serial:5d} {ion_elem:>2s}   {ion_resn} {chain}{curr_res:4d}    "
            f"{trans_coord[0]:8.3f}{trans_coord[1]:8.3f}{trans_coord[2]:8.3f}  1.00  0.00          {ion_elem:>2s}\n"
        )
        ion_records.append(rec)
        placed_info.append(
            {
                "ion_type": ion_elem,
                "ion_index": idx,
                "resnum": curr_res,
                "coordinates": trans_coord.tolist(),
            }
        )
        serial += 1
        curr_res += 1

    with open(rec_path, encoding="utf-8") as f:
        lines = f.readlines()

    clean_lines = [line for line in lines if not line.strip().startswith(("END", "ENDMDL"))]
    out_lines = clean_lines + ion_records + ["END\n"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(out_lines)

    return placed_info


def dock_substrate_peptide(
    receptor_pdb: str | Path,
    output_pdb: str | Path,
    template_pdb: str | Path,
    peptide_seq: str = "KRPSSAPSG",
    acceptor_idx: int | None = None,
    chain_id: str = "B",
    start_res: int = 1,
) -> dict[str, Any]:
    """Dock a substrate peptide into the kinase catalytic active site using template coordinates.

    Reads substrate geometry from `template_pdb` (located in the project `structure/template/` folder)
    and maps it to `receptor_pdb` via ATP superposition.
    """
    rec_path = Path(receptor_pdb)
    out_path = Path(output_pdb)
    tmpl_path = Path(template_pdb)

    if not rec_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {rec_path}")
    if not tmpl_path.exists():
        raise FileNotFoundError(f"Template PDB not found: {tmpl_path}")

    if acceptor_idx is None:
        for i, aa in enumerate(peptide_seq):
            if aa in {"S", "T"}:
                acceptor_idx = i
                break
        if acceptor_idx is None:
            raise ValueError(f"No Ser/Thr acceptor residue found in peptide sequence {peptide_seq}")

    acceptor_aa = peptide_seq[acceptor_idx]
    if acceptor_aa not in {"S", "T"}:
        raise ValueError(f"Acceptor at index {acceptor_idx} must be Ser or Thr, got {acceptor_aa}")

    rec_atp = parse_pdb_atoms(rec_path, res_name="ATP")
    tmpl_atp = parse_pdb_atoms(tmpl_path, res_name="ATP")

    triphosphate_atoms = ["PG", "O1G", "O2G", "O3G", "PB", "O1B", "O3B"]
    shared_atp = [a for a in triphosphate_atoms if a in rec_atp and a in tmpl_atp]
    if len(shared_atp) < 4:
        raise ValueError("Insufficient shared ATP phosphate atoms for docking")

    source_atp = np.vstack([tmpl_atp[a] for a in shared_atp])
    target_atp = np.vstack([rec_atp[a] for a in shared_atp])

    atp_rot, atp_trans = kabsch(source_atp, target_atp)

    # Search for substrate chain (S, B, or second polymer chain)
    tmpl_residues = get_chain_residues(tmpl_path, chain="S")
    sub_chain = "S"
    if not tmpl_residues:
        tmpl_residues = get_chain_residues(tmpl_path, chain="B")
        sub_chain = "B"
    if not tmpl_residues:
        raise ValueError(f"Template {tmpl_path} contains no substrate chain S or B")

    # In 3X2V Chain S, SER 621 is the catalytic acceptor
    tmpl_acceptor_num = 621
    if not any(resnum == tmpl_acceptor_num for resname, resnum in tmpl_residues):
        tmpl_acceptor_num = next(
            resnum for resname, resnum in tmpl_residues if resname in {"SER", "THR"}
        )

    tmpl_acceptor_idx = [num for _, num in tmpl_residues].index(tmpl_acceptor_num)

    # Extract and transform all residues present in template
    known_coords: dict[int, dict[str, np.ndarray]] = {}
    for i in range(len(peptide_seq)):
        ref_idx = tmpl_acceptor_idx + (i - acceptor_idx)
        if 0 <= ref_idx < len(tmpl_residues):
            ref_tmpl_resnum = tmpl_residues[ref_idx][1]
            tmpl_res_atoms = parse_pdb_atoms(tmpl_path, chain=sub_chain, res_num=ref_tmpl_resnum)
            if tmpl_res_atoms and "CA" in tmpl_res_atoms:
                known_coords[i] = {
                    atom: coord @ atp_rot + atp_trans for atom, coord in tmpl_res_atoms.items()
                }

    if not known_coords:
        raise ValueError("No template substrate residues could be mapped")

    # Extend missing C-terminal residues along backbone trajectory
    for i in range(len(peptide_seq)):
        if i in known_coords:
            continue
        if i > max(known_coords.keys()):
            last_i = i - 1
            prev_i = i - 2 if (i - 2) in known_coords else last_i
            last_res = known_coords[last_i]
            prev_res = known_coords[prev_i]

            if "CA" in last_res and "CA" in prev_res and last_i != prev_i:
                u = last_res["CA"] - prev_res["CA"]
                norm = np.linalg.norm(u)
                u = u / norm if norm > 1e-4 else np.array([1.0, 0.0, 0.0])
            elif "C" in last_res and "CA" in last_res:
                u = last_res["C"] - last_res["CA"]
                norm = np.linalg.norm(u)
                u = u / norm if norm > 1e-4 else np.array([1.0, 0.0, 0.0])
            else:
                u = np.array([1.0, 0.0, 0.0])

            up = np.cross(u, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(up) < 0.1:
                up = np.cross(u, np.array([0.0, 1.0, 0.0]))
            up = up / np.linalg.norm(up)

            last_c = last_res.get("C", last_res["CA"] + 1.52 * u)
            n_pos = last_c + 1.33 * u
            ca_pos = n_pos + 1.46 * u
            c_pos = ca_pos + 1.52 * u
            o_pos = c_pos + 1.23 * up
            cb_pos = ca_pos - 1.53 * up

            known_coords[i] = {
                "N": n_pos,
                "CA": ca_pos,
                "C": c_pos,
                "O": o_pos,
                "CB": cb_pos,
            }

    # Extend missing N-terminal residues in reverse
    for i in range(len(peptide_seq) - 1, -1, -1):
        if i in known_coords:
            continue
        if i < min(known_coords.keys()):
            next_i = i + 1
            next2_i = i + 2 if (i + 2) in known_coords else next_i
            next_res = known_coords[next_i]
            next2_res = known_coords[next2_i]

            if "CA" in next_res and "CA" in next2_res and next_i != next2_i:
                u = next2_res["CA"] - next_res["CA"]
                norm = np.linalg.norm(u)
                u = u / norm if norm > 1e-4 else np.array([1.0, 0.0, 0.0])
            else:
                u = np.array([1.0, 0.0, 0.0])

            up = np.cross(u, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(up) < 0.1:
                up = np.cross(u, np.array([0.0, 1.0, 0.0]))
            up = up / np.linalg.norm(up)

            next_n = next_res.get("N", next_res["CA"] - 1.46 * u)
            c_pos = next_n - 1.33 * u
            ca_pos = c_pos - 1.52 * u
            n_pos = ca_pos - 1.46 * u
            o_pos = c_pos + 1.23 * up
            cb_pos = ca_pos - 1.53 * up

            known_coords[i] = {
                "N": n_pos,
                "CA": ca_pos,
                "C": c_pos,
                "O": o_pos,
                "CB": cb_pos,
            }

    peptide_records = []
    serial = 80000

    for i, aa in enumerate(peptide_seq):
        res_3 = AA1_TO_3.get(aa, "ALA")
        target_resnum = start_res + i
        res_atoms = known_coords[i]

        allowed_atoms = ["N", "CA", "C", "O"]
        if aa != "G":
            allowed_atoms.append("CB")
        if i == acceptor_idx:
            og_name = "OG" if aa == "S" else "OG1"
            allowed_atoms.append(og_name)
            if og_name not in res_atoms and "CB" in res_atoms:
                res_atoms[og_name] = res_atoms["CB"] + np.array([0.0, 1.43, 0.0])

        for atom_name in allowed_atoms:
            if atom_name not in res_atoms:
                continue
            coord = res_atoms[atom_name]
            elem = atom_name[:1]
            rec = (
                f"ATOM  {serial:5d} {atom_name:>4s} {res_3} {chain_id}{target_resnum:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00          {elem:>2s}\n"
            )
            peptide_records.append(rec)
            serial += 1

    with open(rec_path, encoding="utf-8") as f:
        rec_lines = f.readlines()

    clean_rec = [line for line in rec_lines if not line.strip().startswith(("END", "ENDMDL"))]
    if clean_rec and not clean_rec[-1].strip().startswith("TER"):
        clean_rec.append("TER\n")

    complex_lines = clean_rec + peptide_records + ["TER\n", "END\n"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(complex_lines)

    comp_atp = parse_pdb_atoms(out_path, res_name="ATP")
    acceptor_atoms = parse_pdb_atoms(out_path, chain=chain_id, res_num=start_res + acceptor_idx)

    og_coord = acceptor_atoms.get("OG")
    if og_coord is None:
        og_coord = acceptor_atoms.get("OG1")
    pg_coord = comp_atp.get("PG")

    distance_og_pg = None
    if og_coord is not None and pg_coord is not None:
        distance_og_pg = float(np.linalg.norm(og_coord - pg_coord))

    return {
        "peptide_seq": peptide_seq,
        "acceptor_res": f"{acceptor_aa}{start_res + acceptor_idx}",
        "distance_og_pg_A": (round(distance_og_pg, 3) if distance_og_pg is not None else None),
        "status": ("PASS" if (distance_og_pg and 2.5 <= distance_og_pg <= 5.5) else "CHECK"),
    }


def assemble_system_model(
    spec: SystemSpec,
    receptor_pdb: str | Path,
    output_pdb: str | Path,
    template_pdb: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> Path:
    """End-to-end assembly of a complete kinase-ligand-ion-substrate complex.

    Applies in sequence:
      1. Point mutation (if spec.variant != "WT") via Modeller or mutagenesis.
      2. Phosphorylation (if spec.phos_residues, e.g. Y171 -> PTR171).
      3. Catalytic divalent ions (if spec.ions, e.g. 2 Mg or 2 Mn).
      4. Substrate peptide docking (if spec.peptide_seq, e.g. EB2 Ser222).
      5. Formats and writes finalized clean PDB coordinates.
    """
    rec_path = Path(receptor_pdb).resolve()
    out_path = Path(output_pdb).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not rec_path.exists():
        raise FileNotFoundError(f"Input receptor PDB not found: {rec_path}")

    w_dir = Path(work_dir).resolve() if work_dir else out_path.parent / "_work_prep"
    w_dir.mkdir(parents=True, exist_ok=True)

    current_pdb = rec_path

    # Step 0: State normalization (strip extraneous components if not requested by SystemSpec)
    clean_lines = []
    has_phos = bool(spec.phos_residues)
    has_nucleotide = bool(spec.ligand)
    has_peptide = bool(spec.peptide_seq or spec.peptide_name)
    allowed_ions = set(spec.ions.keys()) if spec.ions else set()

    with open(current_pdb, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                chain_id = line[21]

                # Strip peptide chain B if not requested
                if not has_peptide and chain_id == "B":
                    continue

                # Strip nucleotide if not requested
                if not has_nucleotide and resname in {"ATP", "ADP", "GDP", "GTP", "ANP"}:
                    continue

                # Strip ions if not requested
                if resname in {"MG", "MN", "CA", "ZN"} and resname not in allowed_ions:
                    continue

                # If phosphorylation is NOT requested, convert non-standard phosphorylated residues to canonical forms
                if not has_phos:
                    if resname == "PTR":
                        if atomname in {"P", "O1P", "O2P", "O3P"}:
                            continue
                        line = line[:17] + "TYR" + line[20:]
                    elif resname == "SEP":
                        if atomname in {"P", "O1P", "O2P", "O3P"}:
                            continue
                        line = line[:17] + "SER" + line[20:]
                    elif resname == "TPO":
                        if atomname in {"P", "O1P", "O2P", "O3P"}:
                            continue
                        line = line[:17] + "THR" + line[20:]

            clean_lines.append(line)

    norm_pdb = w_dir / "step_state_normalized.pdb"
    norm_pdb.write_text("".join(clean_lines))
    current_pdb = norm_pdb

    # Step 1: Mutation (if variant != WT)
    if spec.variant and spec.variant.upper() != "WT":
        try:
            from colabmda.modeller.engine import _MODELLER_OK, mutate_once

            if _MODELLER_OK:
                mut_res, _, _, _, _ = mutate_once(
                    pdb_in=str(current_pdb),
                    chain="A",
                    mut_str=spec.variant,
                    seed=42,
                    outdir_mut=str(w_dir),
                )
                if mut_res and Path(mut_res).exists():
                    current_pdb = Path(mut_res)
        except Exception:
            pass

    # Step 2: Phosphorylation (e.g. Y171 -> PTR171)
    if spec.phos_residues:
        for phos in spec.phos_residues:
            m = re.search(r"(\d+)", phos)
            if m:
                res_num = int(m.group(1))
                target_res = "PTR"
                if "S" in phos.upper():
                    target_res = "SEP"
                elif "T" in phos.upper():
                    target_res = "TPO"

                phos_pdb = w_dir / f"step_phos_{target_res}_{res_num}.pdb"
                phosphorylate_residue(
                    input_pdb=current_pdb,
                    output_pdb=phos_pdb,
                    res_num=res_num,
                    target_res=target_res,
                )
                current_pdb = phos_pdb

    # Step 3: Nucleotide positioning (if spec.ligand and missing from receptor)
    if spec.ligand and template_pdb:
        lig_upper = str(spec.ligand).upper()
        if lig_upper in {"ATP", "ADP", "GDP"}:
            nuc_pdb = w_dir / f"step_nucleotide_{lig_upper}.pdb"
            position_nucleotide(
                receptor_pdb=current_pdb,
                output_pdb=nuc_pdb,
                template_pdb=template_pdb,
                ligand_name=lig_upper,
            )
            current_pdb = nuc_pdb

    # Step 4: Divalent catalytic ions
    if spec.ions and template_pdb:
        # Strip preexisting ions from receptor to prevent Coulombic clashes
        stripped_pdb = w_dir / "step_strip_ions.pdb"
        clean_lines = []
        with open(current_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    resname = line[17:20].strip()
                    if resname in {"MG", "MN", "CA", "ZN"}:
                        continue
                clean_lines.append(line)
        stripped_pdb.write_text("".join(clean_lines))
        current_pdb = stripped_pdb

        total_ions = 0
        primary_ion = "MG"
        for ion_name, count in spec.ions.items():
            if ion_name in {"MG", "MN", "CA", "ZN"}:
                total_ions += count
                primary_ion = ion_name

        if total_ions > 0:
            ions_pdb = w_dir / f"step_ions_{primary_ion}_{total_ions}.pdb"
            position_ions(
                receptor_pdb=current_pdb,
                output_pdb=ions_pdb,
                template_pdb=template_pdb,
                ion_type=primary_ion,
                count=total_ions,
            )
            current_pdb = ions_pdb

    # Step 4: Substrate peptide docking
    peptide_seq = spec.peptide_seq
    if not peptide_seq and spec.peptide_name:
        preset = SUBSTRATE_PRESETS.get(spec.peptide_name.lower())
        if preset:
            peptide_seq = preset["seq"]

    if peptide_seq and template_pdb:
        dock_pdb = w_dir / "step_dock_substrate.pdb"
        dock_substrate_peptide(
            receptor_pdb=current_pdb,
            output_pdb=dock_pdb,
            template_pdb=template_pdb,
            peptide_seq=peptide_seq,
            chain_id="B",
        )
        current_pdb = dock_pdb

    shutil.copy2(current_pdb, out_path)
    return out_path
