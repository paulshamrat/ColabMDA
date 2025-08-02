# OpenMM Protein-Water MD Workflow

This folder provides a streamlined workflow for running molecular dynamics (MD) simulations of proteins in explicit water using OpenMM and PDBFixer. The workflow is designed for ease of use and reproducibility, from PDB download and cleaning to full MD production runs.

## Installation

### 1. Notebook Setup (Colab)
1. Open a new Colab notebook in Google Drive.
2. In the first cell, mount Google Drive and check GPU allocation:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   !nvidia-smi
   ```
   *(This ensures your session has access to Google Drive and a GPU is available.)*

> **Note:** All environment setup and package installation (including Conda, Mamba, OpenMM, and analysis libraries) should be performed in the Colab Terminal, not in notebook cells.

### 2. Installation on Terminal

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
## 00 installation
#Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

#Step 2: Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

#Step 3: Install Mamba into the base environment
conda install -y -n base -c conda-forge mamba

#Step 4: Install CUDA-enabled OpenMM and OpenMMTools
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools

#Step 5: Install PDBFixer (conda, fallback to pip)
conda install -y -c conda-forge pdbfixer || pip install pdbfixer

#Step 6: Install MDAnalysis, MDTraj, NumPy, Matplotlib, and Biopython
mamba install -y -c conda-forge mdanalysis mdtraj numpy matplotlib biopython

#Step 7: Verify installations
python3 - << 'EOF'
from openmm import Platform; print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import MDAnalysis, mdtraj, Bio; print("MDAnalysis:", MDAnalysis.__version__, "MDTraj:", mdtraj.__version__, "Biopython:", Bio.__version__)
EOF
```


## Workflow Overview

1. **Preprocess the PDB Structure**
   - Download and clean a PDB structure, removing heterogens (except water), building missing residues/atoms, and adding hydrogens at pH 7.0.
   - Command:
     ```bash
     python3 pdbfixer_cleaning.py <pdbid>
     # Example:
     python3 pdbfixer_cleaning.py 4ldj
     ```
   - Output: Creates a folder `<pdbid>/` containing `<pdbid>.pdb` (raw) and `<pdbid>_cleaned.pdb` (processed).


2. **Run MD Simulation**
   - Perform minimization, equilibration (NVT & NPT), and production MD in the same folder.
   - Command:
     ```bash
     python3 openmm_proteinwater.py <pdbid> \
         --total-ns 1 \
         --traj-interval 10.0 \
         --equil-time 10.0 \
         --checkpoint-ps 10.0
     # Example:
     python3 openmm_proteinwater.py 4ldj --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
     ```
   - Output: Trajectory files (`.dcd`), logs, checkpoints, and final system files in `<pdbid>/`.


**Note on Colab GPU Allocation:**

Colab sessions may lose GPU allocation after some time, interrupting long MD simulations. If this happens:

1. Re-do the Google Colab setup (mount Google Drive, check GPU).
2. Re-run the installation steps in the Colab Terminal (if needed).
3. Run the same MD command:
   ```bash
     python3 openmm_proteinwater.py <pdbid> \
         --total-ns 1 \
         --traj-interval 10.0 \
         --equil-time 10.0 \
         --checkpoint-ps 10.0
     # Example:
     python3 openmm_proteinwater.py 4ldj --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
   ```
   The simulation will automatically resume from the last checkpoint file (`prod.chk`).

### Adjusting Simulation Parameters

You can increase the simulation length, change how often trajectory frames and checkpoints are written, and adjust equilibration time using the command-line options:

- `--total-ns`: **Total production time in nanoseconds.**
  - Example: `--total-ns 100` for a 100 ns simulation.
- `--traj-interval`: **Trajectory write interval in picoseconds.**
  - Example: `--traj-interval 100` writes a frame every 100 ps.
- `--equil-time`: **Equilibration phase duration in picoseconds.**
  - Example: `--equil-time 100` for longer equilibration.
- `--checkpoint-ps`: **Checkpoint interval in picoseconds.**
  - Example: `--checkpoint-ps 1000` saves a checkpoint every 1000 ps.

**To run a longer simulation (e.g., 100 ns) with less frequent trajectory and checkpoint writes:**
```bash
python3 openmm_proteinwater.py 4ldj \
    --total-ns 100 \
    --traj-interval 100 \
    --equil-time 100 \
    --checkpoint-ps 1000
```
This will run a 100 ns simulation, write trajectory frames every 100 ps, equilibrate for 100 ps per phase, and save checkpoints every 1000 ps.

---
---
For questions or issues, please contact the repository maintainer.
