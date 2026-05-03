#!/usr/bin/env python3
"""
openmm_trajanalysis.py

Analyze trajectory (RMSD, Rg, RMSF) for a single system or multiple replicas.
Automatically detects r1, r2, r3... subfolders if present.

Usage:
    colabmda openmm analysis --name 4ldj_wt
"""

import os, sys, re, glob, argparse, datetime
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import mdtraj as md
import matplotlib.pyplot as plt

def _data_lines(path):
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == '#' or s[0].isalpha(): continue
            yield s

def _infer_interval(sim_dir, n_frames):
    merged_log = os.path.join(sim_dir, "prod_full.log")
    if os.path.isfile(merged_log):
        times = []
        for s in _data_lines(merged_log):
            parts = s.split(',') if ',' in s else s.split()
            if len(parts) >= 2:
                try: times.append(float(parts[1]))
                except: continue
        if len(times) >= 2:
            return (times[-1] - times[0]) / max(1, n_frames - 1)
    return None

def analyze_single(sim_dir, topo_path, traj_path, interval_user=None):
    u = mda.Universe(topo_path, traj_path)
    n_frames = len(u.trajectory)
    interval = interval_user or _infer_interval(sim_dir, n_frames) or 10.0
    # Correct time: Frames are at 1*int, 2*int ... N*int
    times = (np.arange(n_frames) + 1) * interval

    # RMSD
    rmsd_calc = rms.RMSD(u, u, select="backbone", ref_frame=0)
    rmsd_calc.run()
    rmsd = rmsd_calc.results.rmsd[:, 2]

    # Rg
    heavy = u.select_atoms("protein and not name H*")
    rg = []
    for ts in u.trajectory:
        coords = heavy.positions
        cog = coords.mean(axis=0)
        rg.append(np.sqrt(((coords - cog)**2).sum(axis=1).mean()))
    rg = np.array(rg)

    # RMSF
    t = md.load(traj_path, top=topo_path)
    ca_idx = t.topology.select("name CA")
    t.superpose(t, 0, atom_indices=ca_idx)
    rmsf = md.rmsf(t, t, atom_indices=ca_idx) * 10.0
    resids = [t.topology.atom(i).residue.resSeq for i in ca_idx]

    return {"times": times, "rmsd": rmsd, "rg": rg, "rmsf": rmsf, "resids": resids}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid", help="Simulation directory")
    parser.add_argument("-t", "--topology")
    parser.add_argument("-x", "--trajectory")
    parser.add_argument("-i", "--interval", type=float)
    parser.add_argument("-o", "--outdir")
    args = parser.parse_args()

    sim_dir = os.path.abspath(args.pdbid)
    outdir = Path(args.outdir or f"analysis_{os.path.basename(sim_dir)}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Detect replicas
    replica_dirs = sorted(glob.glob(os.path.join(sim_dir, "r[0-9]*")))
    results = []
    labels = []
    
    if not replica_dirs:
        # Single mode
        topo = args.topology or os.path.join(sim_dir, "solvated.pdb")
        traj = args.trajectory or os.path.join(sim_dir, "prod_full.dcd")
        if os.path.exists(traj):
            results.append(analyze_single(sim_dir, topo, traj, args.interval))
            labels.append("Sim")
            # For single mode, save directly in outdir
            res_outdirs = [outdir]
    else:
        # Multi-replica mode
        print(f"Detected {len(replica_dirs)} replicas: {[os.path.basename(r) for r in replica_dirs]}")
        res_outdirs = []
        for rd in replica_dirs:
            topo = os.path.join(rd, "solvated.pdb")
            traj = os.path.join(rd, "prod_full.dcd")
            if os.path.exists(traj):
                results.append(analyze_single(rd, topo, traj, args.interval))
                rep_name = os.path.basename(rd)
                labels.append(rep_name)
                # Individual replica output folder
                rep_out = outdir / rep_name
                rep_out.mkdir(parents=True, exist_ok=True)
                res_outdirs.append(rep_out)
                
                # Copy equilibration QC plot if it exists in the replica folder
                qc_src = Path(rd) / "equilibration_qc.png"
                if qc_src.exists():
                    import shutil
                    shutil.copy2(qc_src, rep_out / "equilibration_qc.png")

    if not results:
        sys.exit("Error: No valid trajectory data found.")

    # Plotting Styles
    plt.rcParams.update({"font.size": 12, "figure.figsize": (10, 6), "axes.grid": True, "grid.alpha": 0.3})
    colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B3', '#CCB974']

    # 1. Individual Plots
    for i, res in enumerate(results):
        target = res_outdirs[i]
        
        # Save CSV Data for Aggregation
        pd.DataFrame({"Time (ps)": res["times"], "RMSD (Å)": res["rmsd"]}).to_csv(target / "rmsd.csv", index=False)
        pd.DataFrame({"Time (ps)": res["times"], "Rg (Å)": res["rg"]}).to_csv(target / "rg.csv", index=False)
        pd.DataFrame({"Residue": res["resids"], "RMSF (Å)": res["rmsf"]}).to_csv(target / "rmsf.csv", index=False)

        # RMSD
        plt.figure(); plt.plot(res["times"], res["rmsd"], color=colors[0])
        plt.xlabel("Time (ps)"); plt.ylabel("RMSD (Å)"); plt.title(f"Backbone RMSD - {labels[i]}")
        plt.savefig(target / "rmsd_vs_time.png", dpi=300); plt.close()
        
        # Rg
        plt.figure(); plt.plot(res["times"], res["rg"], color=colors[1])
        plt.xlabel("Time (ps)"); plt.ylabel("Rg (Å)"); plt.title(f"Radius of Gyration - {labels[i]}")
        plt.savefig(target / "rg_vs_time.png", dpi=300); plt.close()
        
        # RMSF
        plt.figure(); plt.plot(res["resids"], res["rmsf"], color=colors[2])
        plt.xlabel("Residue Index"); plt.ylabel("RMSF (Å)"); plt.title(f"Residue Flexibility - {labels[i]}")
        plt.savefig(target / "rmsf_per_residue.png", dpi=300); plt.close()

    # 2. Aggregated Plots (if multiple replicas)
    if len(results) > 1:
        agg_dir = outdir / "aggregate"
        agg_dir.mkdir(parents=True, exist_ok=True)
        
        # RMSD Avg
        plt.figure()
        min_len = min(len(r["rmsd"]) for r in results)
        rmsd_matrix = np.array([r["rmsd"][:min_len] for r in results])
        avg_rmsd = np.mean(rmsd_matrix, axis=0)
        std_rmsd = np.std(rmsd_matrix, axis=0)
        t_avg = results[0]["times"][:min_len]
        
        # Save Aggregate CSV
        pd.DataFrame({"Time (ps)": t_avg, "Mean": avg_rmsd, "Std": std_rmsd}).to_csv(agg_dir / "rmsd_avg.csv", index=False)

        plt.plot(t_avg, avg_rmsd, color='black', lw=2, label="Average")
        plt.fill_between(t_avg, avg_rmsd - std_rmsd, avg_rmsd + std_rmsd, color='gray', alpha=0.2)
        for i, res in enumerate(results):
            plt.plot(res["times"][:min_len], res["rmsd"][:min_len], alpha=0.3, label=labels[i], color=colors[i%len(colors)])
        plt.xlabel("Time (ps)"); plt.ylabel("RMSD (Å)"); plt.legend(); plt.title("Aggregate Backbone RMSD")
        plt.savefig(agg_dir / "rmsd_aggregate.png", dpi=300); plt.close()

        # Rg Avg
        plt.figure()
        rg_matrix = np.array([r["rg"][:min_len] for r in results])
        avg_rg = np.mean(rg_matrix, axis=0)
        std_rg = np.std(rg_matrix, axis=0)
        
        # Save Aggregate CSV
        pd.DataFrame({"Time (ps)": t_avg, "Mean": avg_rg, "Std": std_rg}).to_csv(agg_dir / "rg_avg.csv", index=False)

        plt.plot(t_avg, avg_rg, color='black', lw=2, label="Average")
        plt.fill_between(t_avg, avg_rg - std_rg, avg_rg + std_rg, color='gray', alpha=0.2)
        for i, res in enumerate(results):
            plt.plot(res["times"][:min_len], res["rg"][:min_len], alpha=0.3, label=labels[i], color=colors[(i+1)%len(colors)])
        plt.xlabel("Time (ps)"); plt.ylabel("Rg (Å)"); plt.legend(); plt.title("Aggregate Radius of Gyration")
        plt.savefig(agg_dir / "rg_aggregate.png", dpi=300); plt.close()

        # RMSF Avg
        plt.figure()
        rmsf_matrix = np.array([r["rmsf"] for r in results])
        avg_rmsf = np.mean(rmsf_matrix, axis=0)
        std_rmsf = np.std(rmsf_matrix, axis=0)
        resids = results[0]["resids"]
        
        # Save Aggregate CSV
        pd.DataFrame({"Residue": resids, "Mean": avg_rmsf, "Std": std_rmsf}).to_csv(agg_dir / "rmsf_avg.csv", index=False)

        plt.plot(resids, avg_rmsf, color='black', lw=2, label="Mean")
        plt.fill_between(resids, avg_rmsf - std_rmsf, avg_rmsf + std_rmsf, alpha=0.2, color='gray', label="Std Dev")
        for i, r in enumerate(rmsf_matrix):
            plt.plot(resids, r, alpha=0.3, label=labels[i], lw=1)
        plt.xlabel("Residue Index"); plt.ylabel("RMSF (Å)"); plt.legend(); plt.title("Aggregate Residue Flexibility")
        plt.savefig(agg_dir / "rmsf_aggregate.png", dpi=300); plt.close()

    print(f"✔ Analysis complete. Nested results and raw data saved in: {outdir}")

if __name__ == "__main__":
    import pandas as pd
    from pathlib import Path
    main()
