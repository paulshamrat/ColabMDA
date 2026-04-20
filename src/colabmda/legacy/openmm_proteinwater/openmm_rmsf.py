#!/usr/bin/env python3
"""
openmm_rmsf.py

Compute C-alpha RMSF and save:
  - rmsf.csv
  - rmsf_per_residue.png
"""

import argparse
import datetime
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md


def parse_args():
    p = argparse.ArgumentParser(description="Compute C-alpha RMSF")
    p.add_argument("simdir", help="Simulation directory (contains solvated.pdb, prod_full.dcd)")
    p.add_argument("-t", "--topology", default=None)
    p.add_argument("-x", "--trajectory", default=None)
    p.add_argument("-o", "--outdir", default=None)
    p.add_argument("--ylim", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    p.add_argument("--resid-min", type=int, default=None, help="Optional lower residue bound")
    p.add_argument("--resid-max", type=int, default=None, help="Optional upper residue bound")
    return p.parse_args()


def main():
    args = parse_args()
    simdir = os.path.abspath(args.simdir)
    label = os.path.basename(simdir.rstrip(os.sep)) or "sim"

    topo = os.path.abspath(args.topology or os.path.join(simdir, "solvated.pdb"))
    traj = os.path.abspath(args.trajectory or os.path.join(simdir, "prod_full.dcd"))
    if not os.path.isfile(topo):
        sys.exit(f"Topology not found: {topo}")
    if not os.path.isfile(traj):
        sys.exit(f"Trajectory not found: {traj}")

    outdir = args.outdir or f"analysis_{label}"
    os.makedirs(outdir, exist_ok=True)

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
        "figure.figsize": (8, 6),
        "savefig.dpi": 600,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
    })
    
    t = md.load(traj, top=topo)
    ca_idx = t.topology.select("name CA")
    t.superpose(t, 0, atom_indices=ca_idx)
    rmsf_A = md.rmsf(t, t, atom_indices=ca_idx) * 10.0
    resids = np.array([t.topology.atom(i).residue.resSeq for i in ca_idx], dtype=int)

    mask = np.ones_like(resids, dtype=bool)
    if args.resid_min is not None:
        mask &= (resids >= args.resid_min)
    if args.resid_max is not None:
        mask &= (resids <= args.resid_max)

    x = resids[mask]
    y = rmsf_A[mask]

    csv_path = os.path.join(outdir, "rmsf.csv")
    np.savetxt(csv_path, np.column_stack([x, y]), delimiter=",", header="residue,rmsf_A", comments="")

    plt.figure()
    plt.plot(x, y, color='firebrick', lw=1.5)
    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title("Residue Flexibility (RMSF)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    if args.ylim is not None:
        plt.ylim(args.ylim[0], args.ylim[1])
    plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"))
    plt.savefig(os.path.join(outdir, "rmsf_per_residue.pdf"))
    plt.close()
    print(f"  → Saved {outdir}/rmsf_per_residue.png (and .pdf)")

    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
