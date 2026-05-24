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

import argparse
import glob
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Merge MD chunk files")
    p.add_argument("pdbid", help="4-letter PDB ID directory")
    p.add_argument(
        "-t", "--topology", default=None, help="Topology PDB (default: <pdbid>/solvated.pdb)"
    )
    p.add_argument("-o", "--out-traj", default="prod_full.dcd", help="Output merged trajectory")
    p.add_argument("-l", "--out-log", default="prod_full.log", help="Output merged log CSV")
    p.add_argument(
        "-s", "--stride", type=int, default=1, help="Keep every Nth frame/row while merging"
    )
    p.add_argument("--center", action="store_true", help="Center the protein in the box")
    p.add_argument(
        "--wrap",
        action="store_true",
        help="Wrap solvent molecules back into the primary box (image_molecules)",
    )
    p.add_argument(
        "--mda", action="store_true", help="Use MDAnalysis for merging and PBC correction"
    )
    p.add_argument("--selection", default=None, help="MDAnalysis selection string (implies --mda)")
    p.add_argument("--ca-only", action="store_true", help="Extract Cα atoms only (implies --mda)")
    p.add_argument(
        "--protein-only", action="store_true", help="Extract protein-only atoms (implies --mda)"
    )
    return p.parse_args()


def extract_start_ps(filename):
    m = re.search(r"prod_(\d+)to\d+ps\.dcd$", filename)
    return int(m.group(1)) if m else float("inf")


def merge_trajectories(pdbid, topology, out_traj, stride=1, center=False, wrap=False):
    import mdtraj as md

    pdbid = os.path.abspath(pdbid)
    out_traj = os.path.abspath(out_traj)
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

    print("Merging DCD chunks (Streaming mode):")
    for f in dcd_files:
        print("  ", f)
    print(f"Using stride: {stride} | center: {center} | wrap: {wrap}")

    # Initialize global frame counter and merged frame count
    global_frame_idx = 0
    merged_count = 0

    # Load topology once to reuse for all chunks
    top_obj = md.load_topology(topo)

    # Open output file for streaming
    total_chunks = len(dcd_files)
    with md.formats.DCDTrajectoryFile(out_traj, "w") as f_out:
        for idx, f in enumerate(dcd_files):
            try:
                # Load one chunk at a time
                chunk = md.load(f, top=top_obj)
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

                    # Write frames to disk immediately
                    f_out.write(
                        sub_chunk.xyz * 10.0,  # convert nm to Angstroms for DCD
                        cell_lengths=(
                            sub_chunk.unitcell_lengths * 10.0
                            if sub_chunk.unitcell_lengths is not None
                            else None
                        ),
                        cell_angles=sub_chunk.unitcell_angles,
                    )
                    merged_count += len(indices)

                global_frame_idx += chunk.n_frames
                print(
                    f" -> Processed chunk {idx + 1}/{total_chunks}: {f} ({chunk.n_frames} frames)"
                )
                sys.stdout.flush()
                # Explicitly delete to free RAM
                del chunk
            except Exception as e:
                print(f"Warning: skipping {f} due to load error: {e}")
                sys.stdout.flush()

    if merged_count == 0:
        os.chdir(start_dir)
        sys.exit("Error: No frames kept after striding. Check your stride value vs total frames.")

    print(f"→ Wrote merged trajectory: {out_traj} ({merged_count} frames)")
    os.chdir(start_dir)


def merge_logs(pdbid, out_log, stride=1):
    pdbid = os.path.abspath(pdbid)
    out_log = os.path.abspath(out_log)
    start_dir = os.getcwd()
    os.chdir(pdbid)
    log_files = sorted(
        glob.glob("prod_*to*ps.log"), key=lambda f: extract_start_ps(f.replace(".log", ".dcd"))
    )
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


def merge_trajectories_mda(
    pdbid,
    topology,
    out_traj,
    stride=1,
    center=False,
    wrap=False,
    selection=None,
    ca_only=False,
    protein_only=False,
):
    try:
        import MDAnalysis as mda_lib
        from MDAnalysis.transformations import center_in_box, unwrap
    except ImportError:
        sys.exit(
            "Error: MDAnalysis is required for --mda options. Install with: pip install mdanalysis"
        )

    pdbid = os.path.abspath(pdbid)
    topo = topology or os.path.join(pdbid, "solvated.pdb")
    topo = os.path.abspath(topo)
    if not os.path.isfile(topo):
        sys.exit(f"Error: topology file not found: {topo}")

    start_dir = os.getcwd()
    os.chdir(pdbid)
    candidates = glob.glob("prod_*to*ps.dcd")
    dcd_files = [f for f in candidates if os.path.isfile(f) and os.path.getsize(f) > 0]
    dcd_files.sort(key=extract_start_ps)
    if not dcd_files:
        os.chdir(start_dir)
        sys.exit("Error: no valid .dcd chunk files found.")

    print("Merging DCD chunks (MDAnalysis mode):")
    for f in dcd_files:
        print("  ", f)
    print(
        f"Using stride: {stride} | center: {center} | wrap: {wrap} | "
        f"ca_only: {ca_only} | protein_only: {protein_only} | selection: {selection}"
    )

    out_traj_abs = os.path.abspath(out_traj)
    out_topo_abs = out_traj_abs.replace(".dcd", ".pdb")

    # Load Universe
    u = mda_lib.Universe(topo, dcd_files)

    # Determine selection
    if ca_only:
        sel_str = "protein and name CA"
    elif protein_only:
        sel_str = "protein"
    elif selection:
        sel_str = selection
    else:
        sel_str = "all"

    output_sel = u.select_atoms(sel_str)
    if len(output_sel) == 0:
        os.chdir(start_dir)
        sys.exit(f"Error: Selection '{sel_str}' matched 0 atoms.")

    # Transformations
    protein = u.select_atoms("protein")
    if (wrap or center) and len(protein) > 0:
        protein.guess_bonds()

    workflow = []
    if wrap and len(protein) > 0:
        workflow.append(unwrap(protein))
    if center and len(protein) > 0:
        workflow.append(center_in_box(protein, center="mass"))

    if workflow:
        u.trajectory.add_transformations(*workflow)

    # Write topology first
    output_sel.write(out_topo_abs)
    print(f"→ Wrote merged topology: {out_topo_abs}")

    # Write trajectory
    total_frames = len(u.trajectory)
    print(f"Writing merged trajectory (total frames: {total_frames})...")
    merged_count = 0
    with mda_lib.Writer(out_traj_abs, output_sel.n_atoms) as W:
        for i, ts in enumerate(u.trajectory):
            if i % stride == 0:
                W.write(output_sel)
                merged_count += 1
            if i % 100 == 0 or i == total_frames - 1:
                print(f" -> Processed frame {i}/{total_frames} ({(i / total_frames) * 100:.1f}%)")
                sys.stdout.flush()

    print(f"→ Wrote merged trajectory: {out_traj_abs} ({merged_count} frames)")
    os.chdir(start_dir)


def main():
    args = parse_args()
    if not os.path.isdir(args.pdbid):
        sys.exit(f"Error: directory '{args.pdbid}' not found.")
    if args.stride < 1:
        sys.exit("Error: --stride must be >= 1")

    # Automatically trigger MDAnalysis mode if any MDAnalysis specific flags are active
    mda_active = args.mda or args.ca_only or args.protein_only or (args.selection is not None)

    if mda_active:
        merge_trajectories_mda(
            args.pdbid,
            args.topology,
            args.out_traj,
            stride=args.stride,
            center=args.center,
            wrap=args.wrap,
            selection=args.selection,
            ca_only=args.ca_only,
            protein_only=args.protein_only,
        )
    else:
        merge_trajectories(
            args.pdbid,
            args.topology,
            args.out_traj,
            stride=args.stride,
            center=args.center,
            wrap=args.wrap,
        )
    merge_logs(args.pdbid, args.out_log, stride=args.stride)


if __name__ == "__main__":
    main()
