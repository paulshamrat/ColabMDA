#!/usr/bin/env python3
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md

# ─── Locate the 1aki subfolder ─────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKDIR  = os.path.join(BASE_DIR, '1aki')

# ─── Paths to logs & trajectory ────────────────────────────────────────────
NVT_LOG   = os.path.join(WORKDIR, 'nvt.log')
NPT_LOG   = os.path.join(WORKDIR, 'npt.log')
PROD_LOG  = os.path.join(WORKDIR, 'prod_1000ps.log')
PDB_FILE  = os.path.join(WORKDIR, 'solvated.pdb')
DCD_FILE  = os.path.join(WORKDIR, 'prod_1000ps.dcd')
OUTDIR    = os.path.join(WORKDIR, 'analysis_plots')
os.makedirs(OUTDIR, exist_ok=True)

# ─── Helper: parse quoted CSV header + data ────────────────────────────────
def read_log(path):
    with open(path, 'r') as f:
        header_line = f.readline().lstrip('#').strip()
    cols = next(csv.reader([header_line]))
    clean = [c.split('(')[0].strip().replace(' ', '_').replace('/', '_') for c in cols]
    df = pd.read_csv(path, skiprows=1, names=clean, sep=',')
    return df

# ─── 1) Equilibration: Temperature & Volume ───────────────────────────────
nvt = read_log(NVT_LOG)   # ['Step','Temperature']
npt = read_log(NPT_LOG)   # ['Step','Temperature','Box_Volume']

for df in (nvt, npt):
    df['Step']    = pd.to_numeric(df['Step'], errors='raise')
    df['Time_ps'] = df['Step'] * 0.002  # 2 fs → 0.002 ps

plt.figure(figsize=(6,4))
plt.plot(nvt['Time_ps'], nvt['Temperature'], label='NVT')
plt.plot(npt['Time_ps'], npt['Temperature'], label='NPT')
plt.xlabel('Time (ps)'); plt.ylabel('Temperature (K)')
plt.title('Temperature vs Time'); plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'temperature.png'))
plt.close()

plt.figure(figsize=(6,4))
plt.plot(npt['Time_ps'], npt['Box_Volume'])
plt.xlabel('Time (ps)'); plt.ylabel('Volume (nm³)')
plt.title('Volume vs Time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'volume.png'))
plt.close()

# ─── 2) Production: Potential Energy ──────────────────────────────────────
prod = read_log(PROD_LOG)
prod['Time']             = pd.to_numeric(prod['Time'], errors='raise')
prod['Potential_Energy'] = pd.to_numeric(prod['Potential_Energy'], errors='raise')
prod['Time_ps']          = prod['Time']

plt.figure(figsize=(6,4))
plt.plot(prod['Time_ps'], prod['Potential_Energy'])
plt.xlabel('Time (ps)'); plt.ylabel('Potential Energy (kJ/mole)')
plt.title('Production Potential Energy vs Time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'potential_energy.png'))
plt.close()

# ─── 3) Trajectory metrics via MDTraj ────────────────────────────────────
traj  = md.load_dcd(DCD_FILE, top=PDB_FILE)
times = traj.time  # ps
bb    = traj.topology.select('backbone')

# 3a) Backbone RMSD vs first frame
rmsd_nm = md.rmsd(traj, traj, frame=0, atom_indices=bb)
plt.figure(figsize=(6,4))
plt.plot(times, rmsd_nm * 10)
plt.xlabel('Time (ps)'); plt.ylabel('RMSD (Å)')
plt.title('Backbone RMSD vs Time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'rmsd.png'))
plt.close()

# 3b) Radius of gyration
rg_nm = md.compute_rg(traj)
plt.figure(figsize=(6,4))
plt.plot(times, rg_nm * 10)
plt.xlabel('Time (ps)'); plt.ylabel('Rg (Å)')
plt.title('Radius of Gyration vs Time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'rg.png'))
plt.close()

# 3c) RMSF per residue
# Use the first frame as reference trajectory for rmsf
reference = traj[0]  # a Trajectory of length 1
rmsf_atom = md.rmsf(traj, atom_indices=bb, reference=reference)
res_map   = [atom.residue.index for atom in traj.topology.atoms if atom.index in bb]
res_ids   = sorted(set(res_map))
rmsf_res  = [
    np.mean([rmsf_atom[i] for i,r in enumerate(res_map) if r == rid]) * 10
    for rid in res_ids
]

plt.figure(figsize=(6,4))
plt.plot(res_ids, rmsf_res)
plt.xlabel('Residue Index'); plt.ylabel('RMSF (Å)')
plt.title('Backbone RMSF per Residue')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'rmsf.png'))
plt.close()

print("✅ Analysis complete! Plots are in:", OUTDIR)
