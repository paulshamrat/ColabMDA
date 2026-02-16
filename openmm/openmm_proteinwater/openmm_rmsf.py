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

    outdir = args.outdir or f"analysis_{label}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(outdir, exist_ok=True)

    t = md.load(traj, top=topo)
    t.superpose(t, 0)
    ca_idx = t.topology.select("name CA")
    rmsf_nm = md.rmsf(t, t, atom_indices=ca_idx) / 10.0
    resids = np.array([t.topology.atom(i).residue.resSeq for i in ca_idx], dtype=int)

    mask = np.ones_like(resids, dtype=bool)
    if args.resid_min is not None:
        mask &= (resids >= args.resid_min)
    if args.resid_max is not None:
        mask &= (resids <= args.resid_max)

    x = resids[mask]
    y = rmsf_nm[mask]

    csv_path = os.path.join(outdir, "rmsf.csv")
    np.savetxt(csv_path, np.column_stack([x, y]), delimiter=",", header="residue,rmsf_nm", comments="")

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, lw=1.5)
    plt.xlabel("Residue")
    plt.ylabel("RMSF (nm)")
    plt.title("Backbone Cα RMSF per Residue")
    plt.grid(True, alpha=0.3)
    if args.ylim is not None:
        plt.ylim(args.ylim[0], args.ylim[1])
    plt.tight_layout()
    fig_path = os.path.join(outdir, "rmsf_per_residue.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()

    print(f"Saved: {csv_path}")
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
