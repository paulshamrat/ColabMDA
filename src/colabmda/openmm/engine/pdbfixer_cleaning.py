#!/usr/bin/env python3
"""
pdbfixer_cleaning.py

Given a PDB ID, this script will:
 1. Create a folder ./<pdbid>/
 2. Download the PDB file into that folder.
 3. Strip out all heterogens except water.
 4. Build any missing residues/atoms.
 5. Add all hydrogens at pH 7.0.
 6. Write out <pdbid>_cleaned.pdb in that folder.

Usage:
    python3 pdbfixer_cleaning.py 4ldj

Requirements:
    conda install -c conda-forge pdbfixer openmm
"""

import argparse
import os
import urllib.request

from openmm.app import PDBFile
from pdbfixer import PDBFixer


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

    from colabmda.openmm.engine.pdbfixer_clean_fromfile import _apply_propka_titration

    success = _apply_propka_titration(input_pdb, temp_fixed_pdb, output_pdb, target_pH, summary_json)

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


def main():
    parser = argparse.ArgumentParser(description="Download & preprocess a PDB by ID")
    parser.add_argument("pdb_id", help="4-character PDB identifier (e.g., 4ldj)")
    parser.add_argument("--outdir", default=None, help="Output directory (default: ./<pdbid>)")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH for titration (default: 7.4)")
    args = parser.parse_args()

    run_clean_by_pdbid(args.pdb_id, args.outdir, ph=args.ph)


if __name__ == "__main__":
    main()
