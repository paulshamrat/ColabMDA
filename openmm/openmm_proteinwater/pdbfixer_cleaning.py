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

import os
import argparse
import urllib.request
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def download_pdb(pdb_id, out_path):
    url = f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb'
    print(f"[Download] Fetching {pdb_id} from RCSB…")
    urllib.request.urlretrieve(url, out_path)
    print(f"[Download] Saved raw PDB → {out_path}")

def preprocess(input_pdb, output_pdb, target_pH=7.0):
    print(f"[Preprocess] Loading {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)
    print("[Preprocess] Stripping heterogens (keeping waters)…")
    fixer.removeHeterogens(keepWater=True)
    print("[Preprocess] Finding/building missing residues & atoms…")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    print(f"[Preprocess] Adding hydrogens at pH {target_pH}…")
    fixer.addMissingHydrogens(pH=target_pH)
    print(f"[Preprocess] Writing cleaned PDB → {output_pdb}")
    with open(output_pdb, 'w') as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

def main():
    parser = argparse.ArgumentParser(description="Download & preprocess a PDB by ID")
    parser.add_argument('pdb_id', help="4-character PDB identifier (e.g., 4ldj)")
    args = parser.parse_args()

    pdb_id = args.pdb_id.lower()
    # Create and enter subdirectory
    os.makedirs(pdb_id, exist_ok=True)
    os.chdir(pdb_id)

    raw_pdb = f"{pdb_id}.pdb"
    cleaned_pdb = f"{pdb_id}_cleaned.pdb"

    # 1) Download raw PDB if missing
    if not os.path.exists(raw_pdb):
        download_pdb(pdb_id, raw_pdb)
    else:
        print(f"[Download] Raw PDB already exists: {raw_pdb}")

    # 2–5) Preprocess
    preprocess(raw_pdb, cleaned_pdb)

    print("✅ All done!")
    print(f"→ Check directory ./{pdb_id}/ for {raw_pdb} and {cleaned_pdb}")

if __name__ == "__main__":
    main()
