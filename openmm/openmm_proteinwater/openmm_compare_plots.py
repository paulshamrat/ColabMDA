#!/usr/bin/env python3
"""
openmm_compare_plots.py

Overlay WT vs mutant metrics from per-system analysis CSV outputs:
  - rmsd.csv (time_ps,rmsd_nm)
  - rg.csv   (time_ps,rg_nm)
  - rmsf.csv (residue,rmsf_nm)

Example:
  python3 openmm_compare_plots.py \
    --series WT=/content/drive/MyDrive/openmm/analysis/single/4ldj_wt \
    --series G12C=/content/drive/MyDrive/openmm/analysis/single/4ldj_G12C \
    --outdir /content/drive/MyDrive/openmm/analysis/compare/4ldj_wt_vs_mutants
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Compare WT and mutant analysis CSVs")
    p.add_argument(
        "--series",
        action="append",
        required=True,
        help="Series mapping in form LABEL=ANALYSIS_DIR (repeatable)",
    )
    p.add_argument("--outdir", required=True, help="Output directory for comparison plots")
    p.add_argument("--rmsd-ylim", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    p.add_argument("--rg-ylim", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    p.add_argument("--rmsf-ylim", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    return p.parse_args()


def parse_series(items):
    pairs = []
    for it in items:
        if "=" not in it:
            raise SystemExit(f"Invalid --series '{it}'. Expected LABEL=DIR")
        label, d = it.split("=", 1)
        label = label.strip()
        d = os.path.abspath(d.strip())
        if not label:
            raise SystemExit(f"Invalid --series '{it}'. Empty label.")
        if not os.path.isdir(d):
            raise SystemExit(f"Analysis dir not found for '{label}': {d}")
        pairs.append((label, d))
    return pairs


def load_csv(path, expected_cols):
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < expected_cols:
        raise SystemExit(f"Unexpected CSV format: {path}")
    return data[:, :expected_cols]


def main():
    args = parse_args()
    series = parse_series(args.series)
    os.makedirs(args.outdir, exist_ok=True)

    # RMSD comparison
    plt.figure(figsize=(10, 6))
    plotted = 0
    for label, d in series:
        arr = load_csv(os.path.join(d, "rmsd.csv"), 2)
        if arr is None:
            continue
        plt.plot(arr[:, 0], arr[:, 1], lw=1.5, label=label)
        plotted += 1
    if plotted:
        plt.xlabel("Time (ps)")
        plt.ylabel("RMSD (nm)")
        plt.title("Backbone RMSD: WT vs Mutants")
        plt.grid(True, alpha=0.3)
        if args.rmsd_ylim is not None:
            plt.ylim(args.rmsd_ylim[0], args.rmsd_ylim[1])
        plt.legend()
        plt.tight_layout()
        out = os.path.join(args.outdir, "compare_rmsd.png")
        plt.savefig(out, dpi=300)
        print(f"Saved: {out}")
    plt.close()

    # Rg comparison
    plt.figure(figsize=(10, 6))
    plotted = 0
    for label, d in series:
        arr = load_csv(os.path.join(d, "rg.csv"), 2)
        if arr is None:
            continue
        plt.plot(arr[:, 0], arr[:, 1], lw=1.5, label=label)
        plotted += 1
    if plotted:
        plt.xlabel("Time (ps)")
        plt.ylabel("Radius of Gyration (nm)")
        plt.title("Rg: WT vs Mutants")
        plt.grid(True, alpha=0.3)
        if args.rg_ylim is not None:
            plt.ylim(args.rg_ylim[0], args.rg_ylim[1])
        plt.legend()
        plt.tight_layout()
        out = os.path.join(args.outdir, "compare_rg.png")
        plt.savefig(out, dpi=300)
        print(f"Saved: {out}")
    plt.close()

    # RMSF comparison
    plt.figure(figsize=(10, 6))
    plotted = 0
    for label, d in series:
        arr = load_csv(os.path.join(d, "rmsf.csv"), 2)
        if arr is None:
            continue
        plt.plot(arr[:, 0], arr[:, 1], lw=1.5, label=label)
        plotted += 1
    if plotted:
        plt.xlabel("Residue")
        plt.ylabel("RMSF (nm)")
        plt.title("Cα RMSF: WT vs Mutants")
        plt.grid(True, alpha=0.3)
        if args.rmsf_ylim is not None:
            plt.ylim(args.rmsf_ylim[0], args.rmsf_ylim[1])
        plt.legend()
        plt.tight_layout()
        out = os.path.join(args.outdir, "compare_rmsf.png")
        plt.savefig(out, dpi=300)
        print(f"Saved: {out}")
    plt.close()

    if not any(os.path.isfile(os.path.join(d, "rmsd.csv")) or
               os.path.isfile(os.path.join(d, "rg.csv")) or
               os.path.isfile(os.path.join(d, "rmsf.csv"))
               for _, d in series):
        sys.exit("No rmsd.csv/rg.csv/rmsf.csv found in provided directories.")


if __name__ == "__main__":
    main()
