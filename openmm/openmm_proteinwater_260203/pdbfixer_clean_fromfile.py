#!/usr/bin/env python3
"""
pdbfixer_clean_fromfile.py

Clean a local PDB file with PDBFixer and write <pdbid>_cleaned.pdb
into an output directory. No downloading.

Usage:
  python3 pdbfixer_clean_fromfile.py --in /path/4LDJ.pdb --outdir /content/work/4ldj_wt --pdbid 4ldj
"""

import os
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB file path (e.g., 4LDJ.pdb)")
    ap.add_argument("--outdir", required=True, help="Output directory (created if missing)")
    ap.add_argument("--pdbid", default="4ldj", help="Folder/name prefix (default: 4ldj)")
    ap.add_argument("--ph", type=float, default=7.0, help="Hydrogen pH (default: 7.0)")
    args = ap.parse_args()

    pdbid = args.pdbid.lower()
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    raw_out = os.path.join(outdir, f"{pdbid}.pdb")
    cleaned_out = os.path.join(outdir, f"{pdbid}_cleaned.pdb")

    # Copy raw PDB into outdir (so everything is self-contained)
    if not os.path.exists(raw_out):
        with open(args.inp, "r") as fin, open(raw_out, "w") as fout:
            fout.write(fin.read())

    print(f"[Preprocess] Loading {raw_out}")
    fixer = PDBFixer(filename=raw_out)
    print("[Preprocess] Removing heterogens (keeping water)...")
    fixer.removeHeterogens(keepWater=True)
    print("[Preprocess] Building missing residues/atoms...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    print(f"[Preprocess] Adding hydrogens at pH {args.ph}...")
    fixer.addMissingHydrogens(pH=args.ph)

    print(f"[Preprocess] Writing cleaned PDB -> {cleaned_out}")
    with open(cleaned_out, "w") as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

    print("✅ Done")
    print("Raw   :", raw_out)
    print("Clean :", cleaned_out)

if __name__ == "__main__":
    main()
