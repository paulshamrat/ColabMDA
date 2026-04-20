#!/usr/bin/env python3
"""
openmm_trajanalysis.py

Compute and plot:
  0) Detected total simulation length
  1) Backbone RMSD vs. time
  2) Radius of gyration vs. time
  3) Cα RMSF per residue

Assumes you have:
  <pdbid>/solvated.pdb
  <pdbid>/prod_full.dcd

Usage:
    python3 openmm_trajanalysis.py <pdbid> [options]

Positional arguments:
  pdbid                    4-letter PDB ID directory (e.g. 4ldj)

Optional arguments:
  -t, --topology PATH      Topology PDB (default: <pdbid>/solvated.pdb)
  -x, --trajectory PATH    Trajectory DCD (default: <pdbid>/prod_full.dcd)
  -i, --interval FLOAT     Frame interval in ps (default: detect from DCD header)
  -o, --outdir DIR         Output directory for plots
"""

import os
import sys
import datetime
import argparse
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import mdtraj as md

def detect_interval_ps(u):
    """Try to read u.trajectory.ts.dt (ps); fallback to None."""
    try:
        return float(u.trajectory.ts.dt)
    except Exception:
        return None

def _data_lines(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and s[0].isdigit():
                yield s

def _infer_interval_from_merged_log(sim_dir, n_frames):
    if n_frames < 2:
        return None
    merged_log = os.path.join(sim_dir, "prod_full.log")
    if not os.path.isfile(merged_log):
        return None
    times = []
    for s in _data_lines(merged_log):
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            times.append(float(parts[1]))
        except Exception:
            continue
    if len(times) < 2:
        return None
    dt = (times[-1] - times[0]) / max(1, len(times) - 1)
    return dt if dt > 0 else None

def _infer_interval_from_chunk_logs(sim_dir):
    logs = glob.glob(os.path.join(sim_dir, "prod_*to*ps.log"))
    if not logs:
        return None
    total_ps = 0.0
    total_frames = 0
    rx = re.compile(r"prod_(\d+)to(\d+)ps\.log$")
    for lp in logs:
        m = rx.search(os.path.basename(lp))
        if not m:
            continue
        start_ps = float(m.group(1))
        end_ps = float(m.group(2))
        total_ps += max(0.0, end_ps - start_ps)
        total_frames += sum(1 for _ in _data_lines(lp))
    if total_ps <= 0 or total_frames <= 0:
        return None
    dt = total_ps / total_frames
    return dt if dt > 0 else None

def parse_args():
    p = argparse.ArgumentParser(description="Analyze merged trajectory for a pdbid")
    p.add_argument("pdbid", help="4-letter PDB ID directory")
    p.add_argument("-t", "--topology",
                   help="Topology PDB (default: <pdbid>/solvated.pdb)")
    p.add_argument("-x", "--trajectory",
                   help="Trajectory DCD (default: <pdbid>/prod_full.dcd)")
    p.add_argument("-i", "--interval", type=float,
                   help="Frame interval in ps (default: detect from header)")
    p.add_argument("-o", "--outdir",
                   help="Directory to save plots (default: analysis_<pdbid>_<timestamp>)")
    return p.parse_args()

def main():
    args = parse_args()
    pdbid = args.pdbid
    sim_dir = os.path.abspath(pdbid)
    sim_label = os.path.basename(sim_dir.rstrip(os.sep)) or "sim"

    # Resolve paths
    top_def   = os.path.join(pdbid, "solvated.pdb")
    traj_def  = os.path.join(pdbid, "prod_full.dcd")
    topo_path = os.path.abspath(args.topology or top_def)
    traj_path = os.path.abspath(args.trajectory or traj_def)

    if not os.path.isfile(topo_path):
        sys.exit(f"Topology file not found: {topo_path}")
    if not os.path.isfile(traj_path):
        sys.exit(f"Trajectory file not found: {traj_path}")

    # Create output directory
    outdir = args.outdir or f"analysis_{sim_label}"
    os.makedirs(outdir, exist_ok=True)

    # Load MDAnalysis Universe
    print("Loading trajectory with MDAnalysis...")
    u = mda.Universe(topo_path, traj_path)
    n_frames = len(u.trajectory)
    # Determine frame spacing
    interval_source = "user"
    interval = args.interval
    if interval is None:
        interval = _infer_interval_from_merged_log(sim_dir, n_frames)
        if interval is not None:
            interval_source = "prod_full.log"
    if interval is None:
        interval = _infer_interval_from_chunk_logs(sim_dir)
        if interval is not None:
            interval_source = "chunk logs/chunk names"
    if interval is None:
        interval = detect_interval_ps(u)
        if interval is not None:
            interval_source = "DCD header"
    if interval is None:
        sys.exit("Failed to detect frame interval; please specify --interval")
    total_ns = (n_frames - 1) * interval / 1000.0
    print(f"Detected {n_frames} frames, interval = {interval:.3f} ps ({interval_source}) → total ≈ {total_ns:.3f} ns\n")

    # Styling for publication
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 10,
        "lines.linewidth": 1.5,
        "figure.figsize": (10, 6),
        "savefig.dpi": 600,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
    })
    
    # 1) RMSD (backbone)
    print("Computing RMSD...")
    rmsd_calc = rms.RMSD(u, u, select="backbone", ref_frame=0)
    rmsd_calc.run()
    rmsd_data = rmsd_calc.results.rmsd  # [frame, time(ps), rmsd(Å), group]
    times     = np.arange(n_frames) * interval
    rmsd_A      = rmsd_data[:,2]

    plt.figure()
    plt.plot(times, rmsd_A, color='#4C72B0', lw=1.0, alpha=1.0)
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (Å)")
    plt.title("Backbone RMSD Stability")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rmsd_vs_time.png"))
    plt.savefig(os.path.join(outdir, "rmsd_vs_time.pdf"))
    plt.close()
    print(f"  → Saved {outdir}/rmsd_vs_time.png (and .pdf)")

    # 2) Radius of gyration
    print("Computing Radius of Gyration...")
    heavy   = u.select_atoms("protein and not name H*")
    rg_vals = []
    for ts in u.trajectory:
        coords = heavy.positions
        cog    = coords.mean(axis=0)
        rg     = np.sqrt(((coords - cog)**2).sum(axis=1).mean())
        rg_vals.append(rg)
    rg_vals = np.array(rg_vals)

    plt.figure()
    plt.plot(times, rg_vals, color='#55A868', lw=1.0, alpha=1.0)
    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Protein Compactness (Rg)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rg_vs_time.png"))
    plt.savefig(os.path.join(outdir, "rg_vs_time.pdf"))
    plt.close()
    print(f"  → Saved {outdir}/rg_vs_time.png (and .pdf)")

    # 3) RMSF (MDTraj)
    print("Computing RMSF...")
    traj = md.load(traj_path, top=topo_path)
    ca_idx  = traj.topology.select("name CA")
    traj.superpose(traj, 0, atom_indices=ca_idx)
    rmsf_A = md.rmsf(traj, traj, atom_indices=ca_idx) * 10.0
    resids  = [traj.topology.atom(i).residue.resSeq for i in ca_idx]

    plt.figure()
    plt.plot(resids, rmsf_A, color='#C44E52', lw=1.5)
    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title("Residue Flexibility (RMSF)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"))
    plt.savefig(os.path.join(outdir, "rmsf_per_residue.pdf"))
    plt.close()
    print(f"  → Saved {outdir}/rmsf_per_residue.png (and .pdf)")

    print(f"\nAll publication-quality plots saved in: {outdir}")

if __name__ == "__main__":
    main()
