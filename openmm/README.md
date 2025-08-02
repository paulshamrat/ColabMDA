
# OpenMM Colab Workflow

This repository explains how to run GPU-accelerated OpenMM simulations in Google Colab, with checkpoint and restart support, and basic trajectory analysis (RMSD, RMSF, Rg).

## Contents

- `00_installation.sh`  
  Step-by-step commands to install Conda, Mamba, OpenMM, PDBFixer, MDAnalysis, MDTraj, Matplotlib, NumPy, and Biopython in the Colab terminal.
- `01_md_openmm.py`  
  Production MD script with automatic checkpointing and restart support.
- `02_trajectory_analysis.py`  
  Analysis script: computes and saves RMSD, per-residue RMSF, and radius of gyration plots.

---

## 1. Notebook Setup

1. **Open a new Colab notebook** in My Drive → New → More → Colab.
2. **In the first cell, mount Google Drive and check GPU:**
    ```python
    from google.colab import drive
    drive.mount('/content/drive')
    !nvidia-smi
    ```
3. **In the second cell, install Conda-for-Colab:**
    ```python
    pip install -q condacolab
    import condacolab
    condacolab.install()
    ```
    *(After this finishes, restart the runtime when prompted.)*

---

## 2. Installation on Terminal

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
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

## 3. Running the Simulation

In the same terminal, run:
```bash
python3 01_openmm_full.py
```
Outputs will be saved to `/content/drive/MyDrive/openmm/1aki/` including:

- `solvated.pdb`
- `system.xml`
- `prod.chk`
- `prod_<ms>ps.dcd`
- `prod_<ms>ps.log`
- `nvt.log`
- `npt.log`

---

## 4. Trajectory Analysis

After frames exist in `prod_<ms>ps.dcd`, run:
```bash
python3 02_trajectory_analysis.py
```
Plots will be saved in `1aki/` as:

- `rmsd.png`
- `rmsf.png`
- `rg.png`

---

## 5. Directory Structure

```text
.
├── 00_installation.sh
├── 01_md_openmm.py
├── 02_trajectory_analysis.py
└── 1aki/
    ├── solvated.pdb
    ├── system.xml
    ├── prod.chk
    ├── prod_<ms>ps.dcd
    ├── prod_<ms>ps.log
    ├── nvt.log
    ├── npt.log
    ├── rmsd.png
    ├── rmsf.png
    └── rg.png
```

---

**Tip:** If you interrupt the simulation, rerun `01_md_openmm.py` to resume from the last checkpoint.
