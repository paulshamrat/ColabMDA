#!/usr/bin/env python3
"""
openmm_rg.py

Compute radius of gyration and save:
  - rg.csv
  - rg_vs_time.png
"""

import argparse
import datetime
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda


def parse_args():
    p = argparse.ArgumentParser(description="Compute radius of gyration")
    p.add_argument("simdir", help="Simulation directory (contains solvated.pdb, prod_full.dcd)")
    p.add_argument("-t", "--topology", default=None)
    p.add_argument("-x", "--trajectory", default=None)
    p.add_argument("-i", "--interval", type=float, default=None, help="Frame interval in ps. Effective = sim_interval * stride")
    p.add_argument("--smooth-ps", type=float, default=25.0, help="Moving-average window in ps")
    p.add_argument("--plot-stride", type=int, default=5, help="Stride for raw line plotting")
    p.add_argument("-o", "--outdir", default=None)
    return p.parse_args()


def _data_lines(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and s[0].isdigit():
                yield s


def detect_interval_from_logs(simdir, n_frames):
    if n_frames > 1:
        merged = os.path.join(simdir, "prod_full.log")
        if os.path.isfile(merged):
            vals = []
            for s in _data_lines(merged):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        vals.append(float(parts[1]))
                    except Exception:
                        pass
            if len(vals) > 1:
                dt = (vals[-1] - vals[0]) / max(1, n_frames - 1)
                if dt > 0:
                    return dt, "prod_full.log"

    logs = glob.glob(os.path.join(simdir, "prod_*to*ps.log"))
    rx = re.compile(r"prod_(\d+)to(\d+)ps\.log$")
    total_ps, total_frames = 0.0, 0
    for lp in logs:
        m = rx.search(os.path.basename(lp))
        if not m:
            continue
        total_ps += max(0.0, float(m.group(2)) - float(m.group(1)))
        total_frames += sum(1 for _ in _data_lines(lp))
    if total_ps > 0 and total_frames > 0:
        dt = total_ps / total_frames
        if dt > 0:
            return dt, "chunk logs/chunk names"
    return None, None


def moving_average(y, win):
    if win <= 1 or win > len(y):
        return None
    kernel = np.ones(win) / win
    return np.convolve(y, kernel, mode="valid")


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

    u = mda.Universe(topo, traj)
    n_frames = len(u.trajectory)

    interval = args.interval
    source = "user"
    if interval is None:
        interval, source = detect_interval_from_logs(simdir, n_frames)
    if interval is None:
        try:
            interval = float(u.trajectory.ts.dt)
            source = "DCD header (WARNING: likely inaccurate)"; print(f"!! WARNING: Using DCD header interval. If total time looks wrong, pass --interval manually.")
        except Exception:
            interval = None
    if interval is None:
        sys.exit("Could not infer interval; pass --interval.")

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
    
    def moving_average(y, win=50):
        if len(y) < win: return y
        return np.convolve(y, np.ones(win)/win, mode='same')

    heavy = u.select_atoms("protein and not name H*")
    rg_vals = []
    for _ts in u.trajectory:
        coords = heavy.positions
        cog = coords.mean(axis=0)
        rg = np.sqrt(((coords - cog) ** 2).sum(axis=1).mean())
        rg_vals.append(rg)
    rg_vals = np.asarray(rg_vals)
    times = np.arange(n_frames) * interval
    rg_smooth = moving_average(rg_vals, win=50)

    csv_path = os.path.join(outdir, "rg.csv")
    np.savetxt(csv_path, np.column_stack([times, rg_vals]), delimiter=",", header="time_ps,rg_A", comments="")

    plt.figure()
    plt.plot(times, rg_vals, color='green', alpha=0.3, label="Raw data")
    plt.plot(times, rg_smooth, color='darkgreen', lw=2, label="Moving average")
    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Protein Compactness (Rg)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rg_vs_time.png"))
    plt.savefig(os.path.join(outdir, "rg_vs_time.pdf"))
    plt.close()
    print(f"  → Saved {outdir}/rg_vs_time.png (and .pdf)")

    print(f"Saved: {csv_path}")
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
