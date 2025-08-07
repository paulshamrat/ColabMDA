#!/usr/bin/env python3
"""
openmm_trajmerge.py

Merge non-overlapping chunk files produced by openmm_proteinwater_chunked* scripts
into one continuous trajectory and one log. Skips any missing or empty chunks.

Usage:
    python3 openmm_trajmerge.py <pdbid> [options]

Positional arguments:
  pdbid                  4-letter PDB ID directory (e.g. 4ldj)

Optional arguments:
  -t, --topology PATH    Path to topology PDB (default: <pdbid>/solvated.pdb)
  -o, --out-traj FILE    Merged trajectory filename (default: prod_full.dcd)
  -l, --out-log FILE     Merged log CSV filename (default: prod_full.log)
"""

import os
import sys
import re
import glob
import argparse

import mdtraj as md
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Merge MD chunk files")
    p.add_argument("pdbid", help="4-letter PDB ID directory")
    p.add_argument("-t", "--topology", default=None,
                   help="Topology PDB (default: <pdbid>/solvated.pdb)")
    p.add_argument("-o", "--out-traj", default="prod_full.dcd",
                   help="Output merged trajectory")
    p.add_argument("-l", "--out-log", default="prod_full.log",
                   help="Output merged log CSV")
    return p.parse_args()

def extract_start_ps(filename):
    m = re.search(r'prod_(\d+)to\d+ps\.dcd$', filename)
    return int(m.group(1)) if m else float('inf')

def merge_trajectories(pdbid, topology, out_traj):
    # determine topology path
    topo = topology or os.path.join(pdbid, "solvated.pdb")
    topo = os.path.abspath(topo)
    if not os.path.isfile(topo):
        sys.exit(f"Error: topology file not found: {topo}")

    os.chdir(pdbid)
    # find candidate chunks
    candidates = glob.glob("prod_*to*ps.dcd")
    # filter to only existing non-empty files
    dcd_files = [f for f in candidates if os.path.isfile(f) and os.path.getsize(f) > 0]
    dcd_files.sort(key=extract_start_ps)
    if not dcd_files:
        sys.exit("Error: no valid .dcd chunk files found.")

    print("Merging DCD chunks:")
    for f in dcd_files:
        print("  ", f)

    # load first
    try:
        merged = md.load(dcd_files[0], top=topo)
    except Exception as e:
        sys.exit(f"Failed to load {dcd_files[0]}: {e}")

    # append the rest
    for f in dcd_files[1:]:
        try:
            chunk = md.load(f, top=topo)
            merged = merged.join(chunk)
        except Exception as e:
            print(f"Warning: skipping {f} due to load error: {e}")

    merged.save_dcd(out_traj)
    print(f"→ Wrote merged trajectory: {out_traj} ({merged.n_frames} frames)")

    os.chdir("..")

def merge_logs(pdbid, out_log):
    os.chdir(pdbid)
    log_files = sorted(glob.glob("prod_*to*ps.log"),
                       key=lambda f: extract_start_ps(f.replace(".log", ".dcd")))
    if not log_files:
        sys.exit("Error: no .log chunk files found.")

    print("Merging log chunks:")
    header_written = False
    with open(out_log, "w") as fout:
        for f in log_files:
            print("  ", f)
            with open(f) as fin:
                for i, line in enumerate(fin):
                    if i == 0:
                        if header_written:
                            continue
                        header_written = True
                    fout.write(line)
    print(f"→ Wrote merged log: {out_log} ({sum(1 for _ in open(out_log))} lines)")

    os.chdir("..")

def main():
    args = parse_args()
    if not os.path.isdir(args.pdbid):
        sys.exit(f"Error: directory '{args.pdbid}' not found.")
    merge_trajectories(args.pdbid, args.topology, args.out_traj)
    merge_logs(args.pdbid, args.out_log)

if __name__ == "__main__":
    main()
