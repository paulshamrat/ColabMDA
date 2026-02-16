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
    p.add_argument("-i", "--interval", type=float, default=None, help="Frame interval in ps")
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
                dt = (vals[-1] - vals[0]) / max(1, len(vals) - 1)
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

    outdir = args.outdir or f"analysis_{label}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
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
            source = "DCD header"
        except Exception:
            interval = None
    if interval is None:
        sys.exit("Could not infer interval; pass --interval.")

    print(f"Frames={n_frames}, interval={interval:.3f} ps ({source})")

    heavy = u.select_atoms("not name H*")
    rg_vals = []
    for _ts in u.trajectory:
        coords = heavy.positions
        cog = coords.mean(axis=0)
        rg = np.sqrt(((coords - cog) ** 2).sum(axis=1).mean()) / 10.0
        rg_vals.append(rg)
    rg_vals = np.asarray(rg_vals)
    times = np.arange(n_frames) * interval

    csv_path = os.path.join(outdir, "rg.csv")
    np.savetxt(csv_path, np.column_stack([times, rg_vals]), delimiter=",", header="time_ps,rg_nm", comments="")

    stride = max(1, args.plot_stride)
    plt.figure(figsize=(10, 6))
    plt.plot(times[::stride], rg_vals[::stride], lw=0.9, alpha=0.4, label=f"Raw (stride={stride})")

    win = max(1, int(round(args.smooth_ps / interval)))
    y_s = moving_average(rg_vals, win)
    if y_s is not None:
        t_s = times[win - 1 :]
        plt.plot(t_s, y_s, lw=2.0, label=f"Moving avg ({args.smooth_ps:g} ps)")

    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (nm)")
    plt.title("Rg vs. Time")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    fig_path = os.path.join(outdir, "rg_vs_time.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()

    print(f"Saved: {csv_path}")
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
