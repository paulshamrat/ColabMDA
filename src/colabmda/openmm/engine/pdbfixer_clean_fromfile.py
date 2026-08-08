#!/usr/bin/env python3
"""
pdbfixer_clean_fromfile.py

Clean a local PDB file with PDBFixer and write <pdbid>_cleaned.pdb
into an output directory. No downloading.

Usage:
  python3 pdbfixer_clean_fromfile.py --in /path/4LDJ.pdb --outdir /content/work/4ldj_wt --pdbid 4ldj
"""

import argparse
import json
import os
import shutil
import subprocess

from openmm.app import PDBFile
from pdbfixer import PDBFixer


def _apply_propka_titration(raw_pdb: str, temp_fixed_pdb: str, out_pdb: str, ph: float, summary_json: str) -> bool:
    """Run PDB2PQR with PROPKA at specified pH to assign explicit AMBER titration states (HID/HIE/HIP/ASH/GLH/LYN/CYX)."""
    pdb2pqr_bin = shutil.which("pdb2pqr")
    if not pdb2pqr_bin:
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
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if res.returncode == 0 and os.path.exists(temp_prot_pdb) and os.path.getsize(temp_prot_pdb) > 0:
            # Audit assigned titratable residue forms
            assigned_forms = {}
            with open(temp_prot_pdb, "r") as fin, open(out_pdb, "w") as fout:
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

            print(f"[Protonation] Successfully assigned {len(assigned_forms)} explicit titratable residue states at pH {ph}.")
            return True
    except Exception as exc:
        print(f"[Protonation] Warning: PDB2PQR titration failed ({exc}). Falling back to PDBFixer heuristics.")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB file path (e.g., 4LDJ.pdb)")
    ap.add_argument("--outdir", required=True, help="Output directory (created if missing)")
    ap.add_argument("--pdbid", default="4ldj", help="Folder/name prefix (default: 4ldj)")
    ap.add_argument("--ph", type=float, default=7.4, help="Hydrogen pH (default: 7.4)")
    args = ap.parse_args()

    run_clean_from_file(args.inp, args.outdir, args.pdbid, args.ph)


if __name__ == "__main__":
    main()
