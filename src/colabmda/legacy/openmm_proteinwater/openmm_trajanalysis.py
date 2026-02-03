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
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or f"analysis_{pdbid}_{timestamp}"
    os.makedirs(outdir, exist_ok=True)

    # Load MDAnalysis Universe
    print("Loading trajectory with MDAnalysis...")
    u = mda.Universe(topo_path, traj_path)
    n_frames = len(u.trajectory)
    # Determine frame spacing
    interval = args.interval or detect_interval_ps(u)
    if interval is None:
        sys.exit("Failed to detect frame interval; please specify --interval")
    total_ns = (n_frames - 1) * interval / 1000.0
    print(f"Detected {n_frames} frames, interval = {interval:.2f} ps → total ≈ {total_ns:.3f} ns\n")

    # 1) RMSD (backbone)
    print("Computing RMSD...")
    rmsd_calc = rms.RMSD(u, u, select="backbone", ref_frame=0)
    rmsd_calc.run()
    rmsd_data = rmsd_calc.rmsd  # [frame, time(ps), rmsd(Å), group]
    times     = np.arange(n_frames) * interval
    rmsd_nm   = rmsd_data[:,2] / 10.0

    plt.figure()
    plt.plot(times, rmsd_nm, "-o", markersize=3)
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (nm)")
    plt.title("Backbone RMSD vs. Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rmsd_vs_time.png"), dpi=300)
    plt.close()
    print(f"  → Saved {outdir}/rmsd_vs_time.png")

    # 2) Radius of gyration
    print("Computing Radius of Gyration...")
    heavy   = u.select_atoms("not name H*")
    rg_vals = []
    for ts in u.trajectory:
        coords = heavy.positions
        cog    = coords.mean(axis=0)
        rg     = np.sqrt(((coords - cog)**2).sum(axis=1).mean()) / 10.0
        rg_vals.append(rg)
    rg_vals = np.array(rg_vals)

    plt.figure()
    plt.plot(times, rg_vals, "-o", markersize=3)
    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (nm)")
    plt.title("Rg vs. Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rg_vs_time.png"), dpi=300)
    plt.close()
    print(f"  → Saved {outdir}/rg_vs_time.png")

    # 3) RMSF (MDTraj)
    print("Computing RMSF...")
    traj = md.load(traj_path, top=topo_path)
    traj.superpose(traj, 0)
    ca_idx  = traj.topology.select("name CA")
    rmsf_A  = md.rmsf(traj, traj, atom_indices=ca_idx)
    rmsf_nm = rmsf_A / 10.0
    resids  = [traj.topology.atom(i).residue.resSeq for i in ca_idx]

    plt.figure()
    plt.plot(resids, rmsf_nm, "-o", markersize=3)
    plt.xlabel("Residue")
    plt.ylabel("RMSF (nm)")
    plt.title("Backbone Cα RMSF per Residue")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"), dpi=300)
    plt.close()
    print(f"  → Saved {outdir}/rmsf_per_residue.png")

    print(f"\nAll plots saved in: {outdir}")

if __name__ == "__main__":
    main()
