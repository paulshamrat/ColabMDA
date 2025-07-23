#!/usr/bin/env python3

import os
# Switch into the folder where your simulation outputs live
os.chdir('1aki')

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# ----------------------
# File names (relative to 1aki/)
pdb_file = 'solvated.pdb'
dcd_file = 'prod_1000ps.dcd'
# ----------------------

# Load whatever frames you have so far
traj = md.load_dcd(dcd_file, top=pdb_file)

# 1) Superpose all frames to the first (backbone only)
backbone = traj.topology.select('backbone')
traj.superpose(traj, 0, atom_indices=backbone)

# 2) Compute RMSD (backbone vs frame 0), in Å
rmsd_nm = md.rmsd(traj, traj, 0, atom_indices=backbone)
rmsd = rmsd_nm * 10.0   # convert nm → Å
times = traj.time       # in ps

# 3) Compute RMSF (per-atom → per-residue), in Å
coords = traj.xyz[:, backbone, :]               # (n_frames, n_atoms, 3)
mean_coords = coords.mean(axis=0)[None, :, :]
flucts = np.sqrt(((coords - mean_coords)**2).mean(axis=0))  # nm per atom

# Map each backbone atom to its residue index
residue_ids = np.array([traj.topology.atom(idx).residue.index for idx in backbone])
unique_res = np.unique(residue_ids)
# Average atom fluctuations to per-residue
rmsf_res = [flucts[residue_ids == rid].mean() * 10.0 for rid in unique_res]  # Å

# 4) Compute radius of gyration (Rg) in nm
rg = md.compute_rg(traj)

# --- Plotting & saving ---

# RMSD vs Time
plt.figure()
plt.plot(times / 1000.0, rmsd)
plt.xlabel("Time (ns)")
plt.ylabel("Backbone RMSD (Å)")
plt.title("RMSD vs Time")
plt.tight_layout()
plt.savefig('rmsd.png')
plt.close()

# Per-Residue RMSF
plt.figure()
plt.plot(unique_res + 1, rmsf_res)  # +1 for 1-based residue numbering
plt.xlabel("Residue Number")
plt.ylabel("RMSF (Å)")
plt.title("Per-Residue RMSF")
plt.tight_layout()
plt.savefig('rmsf.png')
plt.close()

# Radius of Gyration vs Time
plt.figure()
plt.plot(times / 1000.0, rg)
plt.xlabel("Time (ns)")
plt.ylabel("Radius of Gyration (nm)")
plt.title("Rg vs Time")
plt.tight_layout()
plt.savefig('rg.png')
plt.close()

print("Analysis complete. Plots saved as rmsd.png, rmsf.png, rg.png in 1aki/") 
