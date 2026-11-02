#!/usr/bin/env python3
"""
prepare.py

Clean a local PDB file with PDBFixer and write <pdbid>_cleaned.pdb
into an output directory. No downloading.

Usage:
  python3 prepare.py --in /path/4LDJ.pdb --outdir /content/work/4ldj_wt --pdbid 4ldj
"""

import json
import os
import shutil
import subprocess
import urllib.request

from openmm.app import PDBFile
from pdbfixer import PDBFixer


def _apply_propka_titration(
    raw_pdb: str, temp_fixed_pdb: str, out_pdb: str, ph: float, summary_json: str
) -> bool:
    """Run PDB2PQR with PROPKA at specified pH to assign explicit AMBER titration states (HID/HIE/HIP/ASH/GLH/LYN/CYX)."""
    pdb2pqr_bin = shutil.which("pdb2pqr")
    if not pdb2pqr_bin:
        for cand in [
            "/root/miniforge3/envs/openmm_env/bin/pdb2pqr",
            os.path.expanduser("~/miniforge3/envs/openmm_env/bin/pdb2pqr"),
        ]:
            if os.path.exists(cand):
                pdb2pqr_bin = cand
                break

    if not pdb2pqr_bin:
        print("[Protonation] Notice: 'pdb2pqr' is not found in the active environment.")
        print(
            "[Protonation] To enable 3D PROPKA pKa titration, run: mamba install -y -c conda-forge pdb2pqr propka"
        )
        return False

    temp_pqr = out_pdb + ".tmp.pqr"
    temp_prot_pdb = out_pdb + ".tmp_prot.pdb"
    cmd = [
        pdb2pqr_bin,
        "--ff=AMBER",
        "--ffout=AMBER",
        "--titration-state-method=propka",
        f"--with-ph={ph}",
        "--drop-water",
        "--pdb-output",
        temp_prot_pdb,
        temp_fixed_pdb,
        temp_pqr,
    ]
    try:
        print(f"[Protonation] Running PDB2PQR with PROPKA titration at pH {ph}...")
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if (
            res.returncode == 0
            and os.path.exists(temp_prot_pdb)
            and os.path.getsize(temp_prot_pdb) > 0
        ):
            # Audit assigned titratable residue forms
            assigned_forms = {}
            with open(temp_prot_pdb) as fin, open(out_pdb, "w") as fout:
                for line in fin:
                    fout.write(line)
                    if line.startswith(("ATOM", "HETATM")):
                        resname = line[17:20].strip()
                        resnum = line[22:26].strip()
                        chain = line[21].strip()
                        key = f"{resname}_{chain}_{resnum}"
                        if resname in {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX"}:
                            assigned_forms[key] = resname

            # Cleanup temp files
            for p in (temp_pqr, temp_prot_pdb):
                if os.path.exists(p):
                    os.remove(p)

            audit_data = {
                "method": "PDB2PQR / PROPKA",
                "pH": ph,
                "explicit_titration_forms": assigned_forms,
            }
            with open(summary_json, "w") as fj:
                json.dump(audit_data, fj, indent=2)

            print(
                f"[Protonation] Successfully assigned {len(assigned_forms)} explicit titratable residue states at pH {ph}."
            )
            return True
    except Exception as exc:
        print(
            f"[Protonation] Warning: PDB2PQR titration failed ({exc}). Falling back to PDBFixer heuristics."
        )
        for p in (temp_pqr, temp_prot_pdb):
            if os.path.exists(p):
                os.remove(p)
    return False


def run_clean_from_file(pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0):
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    raw_out = os.path.join(outdir, f"{pdbid}.pdb")
    cleaned_out = os.path.join(outdir, f"{pdbid}_cleaned.pdb")
    summary_json = os.path.join(outdir, "protonation_summary.json")

    # Copy raw PDB into outdir (so everything is self-contained)
    if not os.path.exists(raw_out):
        with open(pdb_file) as fin, open(raw_out, "w") as fout:
            fout.write(fin.read())

    print(f"[Preprocess] Loading {raw_out}")
    fixer = PDBFixer(filename=raw_out)
    print("[Preprocess] Removing heterogens (keeping water)...")
    fixer.removeHeterogens(keepWater=True)
    print("[Preprocess] Building missing residues/atoms...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Intermediate fixed PDB without hydrogens
    temp_fixed_pdb = cleaned_out + ".tmp_noH.pdb"
    with open(temp_fixed_pdb, "w") as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

    # Attempt PDB2PQR + PROPKA titration first
    success = _apply_propka_titration(raw_out, temp_fixed_pdb, cleaned_out, ph, summary_json)

    if not success:
        print(f"[Preprocess] Adding hydrogens using PDBFixer heuristics at pH {ph}...")
        fixer.addMissingHydrogens(pH=ph)
        with open(cleaned_out, "w") as out:
            PDBFile.writeFile(fixer.topology, fixer.positions, out)
        with open(summary_json, "w") as fj:
            json.dump({"method": f"PDBFixer heuristics (pH {ph})", "pH": ph}, fj, indent=2)

    if os.path.exists(temp_fixed_pdb):
        os.remove(temp_fixed_pdb)

    print("✅ Done")
    print("Raw    :", raw_out)
    print("Clean  :", cleaned_out)
    print("Summary:", summary_json)


#!/usr/bin/env python3
"""
prepare.py

Given a PDB ID, this script will:
 1. Create a folder ./<pdbid>/
 2. Download the PDB file into that folder.
 3. Strip out all heterogens except water.
 4. Build any missing residues/atoms.
 5. Add all hydrogens at pH 7.0.
 6. Write out <pdbid>_cleaned.pdb in that folder.

Usage:
    python3 prepare.py 4ldj

Requirements:
    conda install -c conda-forge pdbfixer openmm
"""


def download_pdb(pdb_id, out_path):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"[Download] Fetching {pdb_id} from RCSB…")
    urllib.request.urlretrieve(url, out_path)
    print(f"[Download] Saved raw PDB → {out_path}")


def preprocess(input_pdb, output_pdb, target_pH=7.4):
    print(f"[Preprocess] Loading {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)
    print("[Preprocess] Stripping heterogens (keeping waters)…")
    fixer.removeHeterogens(keepWater=True)
    print("[Preprocess] Finding/building missing residues & atoms…")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    outdir = os.path.dirname(os.path.abspath(output_pdb))
    summary_json = os.path.join(outdir, "protonation_summary.json")

    # Intermediate fixed PDB without hydrogens
    temp_fixed_pdb = output_pdb + ".tmp_noH.pdb"
    with open(temp_fixed_pdb, "w") as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

    success = _apply_propka_titration(
        input_pdb, temp_fixed_pdb, output_pdb, target_pH, summary_json
    )

    if not success:
        print(f"[Preprocess] Adding hydrogens using PDBFixer heuristics at pH {target_pH}…")
        fixer.addMissingHydrogens(pH=target_pH)
        print(f"[Preprocess] Writing cleaned PDB → {output_pdb}")
        with open(output_pdb, "w") as out:
            PDBFile.writeFile(fixer.topology, fixer.positions, out)

    if os.path.exists(temp_fixed_pdb):
        os.remove(temp_fixed_pdb)


def run_clean_by_pdbid(pdb_id: str, outdir: str | None = None, ph: float = 7.4):
    pdb_id = pdb_id.lower()
    target_dir = os.path.abspath(outdir or pdb_id)
    os.makedirs(target_dir, exist_ok=True)

    start_dir = os.getcwd()
    try:
        os.chdir(target_dir)
        raw_pdb = f"{pdb_id}.pdb"
        cleaned_pdb = f"{pdb_id}_cleaned.pdb"

        # 1) Download raw PDB if missing
        if not os.path.exists(raw_pdb):
            download_pdb(pdb_id, raw_pdb)
        else:
            print(f"[Download] Raw PDB already exists: {raw_pdb}")

        # 2–5) Preprocess with target pH
        preprocess(raw_pdb, cleaned_pdb, target_pH=ph)

        print("✅ All done!")
        print(f"→ Check directory {target_dir} for {raw_pdb} and {cleaned_pdb}")
    finally:
        os.chdir(start_dir)
