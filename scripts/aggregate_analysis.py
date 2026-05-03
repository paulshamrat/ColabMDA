#!/usr/bin/env python3
import os, sys, argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_data(csv_path):
    if not os.path.exists(csv_path): return None
    return pd.read_csv(csv_path)

def aggregate_system(system_dirs, metric_file):
    dfs = []
    for d in system_dirs:
        df = load_data(os.path.join(d, metric_file))
        if df is not None: dfs.append(df)
    
    if not dfs: return None, None, None
    
    # Align by index/time
    combined = pd.concat(dfs)
    col_name = dfs[0].columns[1] # e.g. 'RMSD (A)' or 'Rg (A)' or 'RMSF (A)'
    time_col = dfs[0].columns[0]
    
    grouped = combined.groupby(combined.iloc[:, 0])
    mean = grouped[col_name].mean()
    std = grouped[col_name].std()
    return mean.index, mean.values, std.values

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--series", action="append", help="LABEL=DIR1,DIR2,DIR3")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Styling
    plt.style.use('seaborn-v0_8-muted')
    plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.3})
    
    metrics = {
        "rmsd.csv": ("Time (ps)", "RMSD (Å)", "System Stability (RMSD)"),
        "rg.csv": ("Time (ps)", "Radius of Gyration (Å)", "Compactness (Rg)"),
        "rmsf.csv": ("Residue Index", "RMSF (Å)", "Flexibility (RMSF)")
    }
    
    for filename, (xlabel, ylabel, title) in metrics.items():
        plt.figure(figsize=(10, 6))
        for item in args.series:
            label, dirs_str = item.split("=")
            dirs = [d.strip() for d in dirs_str.split(",")]
            
            x, mean, std = aggregate_system(dirs, filename)
            if x is None: continue
            
            p = plt.plot(x, mean, label=f"{label} (avg)", lw=2)
            color = p[0].get_color()
            plt.fill_between(x, mean-std, mean+std, color=color, alpha=0.2)
            
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, filename.replace(".csv", "_avg.png")), dpi=300)
        plt.close()
    
    print(f"✅ Aggregate plots saved in {args.outdir}")

if __name__ == "__main__":
    main()
