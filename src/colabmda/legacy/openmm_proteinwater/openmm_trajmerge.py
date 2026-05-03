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
  -s, --stride INT       Keep every Nth frame/row while merging (default: 1)
"""

import os
import sys
import re
import glob
import argparse

import mdtraj as md

def parse_args():
    p = argparse.ArgumentParser(description="Merge MD chunk files")
    p.add_argument("pdbid", help="4-letter PDB ID directory")
    p.add_argument("-t", "--topology", default=None,
                   help="Topology PDB (default: <pdbid>/solvated.pdb)")
    p.add_argument("-o", "--out-traj", default="prod_full.dcd",
                   help="Output merged trajectory")
    p.add_argument("-l", "--out-log", default="prod_full.log",
                   help="Output merged log CSV")
    p.add_argument("-s", "--stride", type=int, default=1,
                   help="Keep every Nth frame/row while merging")
    p.add_argument("--center", action="store_true",
                   help="Center the protein in the box")
    p.add_argument("--wrap", action="store_true",
                   help="Wrap solvent molecules back into the primary box (image_molecules)")
    return p.parse_args()

def extract_start_ps(filename):
    m = re.search(r'prod_(\d+)to\d+ps\.dcd$', filename)
    return int(m.group(1)) if m else float('inf')

def merge_trajectories(pdbid, topology, out_traj, stride=1, center=False, wrap=False):
    pdbid = os.path.abspath(pdbid)
    # determine topology path
    topo = topology or os.path.join(pdbid, "solvated.pdb")
    topo = os.path.abspath(topo)
    if not os.path.isfile(topo):
        sys.exit(f"Error: topology file not found: {topo}")

    start_dir = os.getcwd()
    os.chdir(pdbid)
    # find candidate chunks
    candidates = glob.glob("prod_*to*ps.dcd")
    # filter to only existing non-empty files
    dcd_files = [f for f in candidates if os.path.isfile(f) and os.path.getsize(f) > 0]
    dcd_files.sort(key=extract_start_ps)
    if not dcd_files:
        os.chdir(start_dir)
        sys.exit("Error: no valid .dcd chunk files found.")

    print("Merging DCD chunks:")
    for f in dcd_files:
        print("  ", f)
    print(f"Using stride: {stride} | center: {center} | wrap: {wrap}")

    # Initialize global frame counter and collection list
    global_frame_idx = 0
    to_join = []

    for f in dcd_files:
        try:
            chunk = md.load(f, top=topo)
            # Find which frames in this chunk belong in the thinned trajectory
            indices = [i for i in range(chunk.n_frames) if (global_frame_idx + i) % stride == 0]
            if indices:
                sub_chunk = chunk[indices]
                if wrap:
                    sub_chunk.image_molecules(inplace=True)
                if center:
                    # Select protein for centering
                    protein_indices = sub_chunk.topology.select("protein")
                    if len(protein_indices) > 0:
                        sub_chunk.center_coordinates()
                to_join.append(sub_chunk)
            global_frame_idx += chunk.n_frames
        except Exception as e:
            print(f"Warning: skipping {f} due to load error: {e}")

    if not to_join:
        os.chdir(start_dir)
        sys.exit("Error: No frames kept after striding. Check your stride value vs total frames.")

    # Merge the selected frames
    merged = to_join[0]
    for subset in to_join[1:]:
        merged = merged.join(subset)

    merged.save_dcd(out_traj)
    print(f"→ Wrote merged trajectory: {out_traj} ({merged.n_frames} frames)")
    os.chdir(start_dir)

def merge_logs(pdbid, out_log, stride=1):
    pdbid = os.path.abspath(pdbid)
    start_dir = os.getcwd()
    os.chdir(pdbid)
    log_files = sorted(glob.glob("prod_*to*ps.log"),
                       key=lambda f: extract_start_ps(f.replace(".log", ".dcd")))
    if not log_files:
        os.chdir(start_dir)
        sys.exit("Error: no .log chunk files found.")

    print("Merging log chunks:")
    print(f"Using stride: {stride}")
    header_written = False
    frame_idx = 0
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
                        continue

                    # Keep every Nth data row to match merged trajectory stride.
                    if frame_idx % max(1, stride) == 0:
                        fout.write(line)
                    frame_idx += 1
    print(f"→ Wrote merged log: {out_log} ({sum(1 for _ in open(out_log))} lines)")
    os.chdir(start_dir)

def main():
    args = parse_args()
    if not os.path.isdir(args.pdbid):
        sys.exit(f"Error: directory '{args.pdbid}' not found.")
    if args.stride < 1:
        sys.exit("Error: --stride must be >= 1")
    merge_trajectories(args.pdbid, args.topology, args.out_traj, stride=args.stride, center=args.center, wrap=args.wrap)
    merge_logs(args.pdbid, args.out_log, stride=args.stride)

if __name__ == "__main__":
    main()
