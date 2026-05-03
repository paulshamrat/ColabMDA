# ColabMDA 🧬
**Publication-grade Molecular Dynamics on Google Colab**

ColabMDA is a modular, high-performance pipeline for protein modeling and MD simulation. It is specifically optimized for Google Colab and Google Drive, with built-in protection against GPU timeouts.

---

## 1. Installation

### 1.1. Mounting Drive (Required)
First, mount your Google Drive to ensure your data is persistent.
```python
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

### 1.2. Environment Setup (Required)
Install OpenMM, Modeller, and GPU drivers. Run this in the **Colab Terminal**.
```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest
```

### 1.3. Package Installation (Required)
Install the ColabMDA command-line tool.
```bash
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Steps 1.2 and 1.3** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your `.chk` files and resume from where it left off.

---

<details>
<summary>Manual / Advanced Installation (Not for general users)</summary>

### Manual Terminal Installation
```bash
# Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o /tmp/miniforge.sh && \
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

# Step 2: Initialize Conda
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Step 3: Install Core Stack
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools mdanalysis mdtraj numpy matplotlib biopython
```
</details>

---

## 🛠 2. Simulation Workflow

### 2.1. Build Structures (WT and Mutants)
**Environment:** `modeller_env` | **Directory:** `/content/drive/MyDrive/openmm`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env
cd /content/drive/MyDrive/openmm

# Example: Build Wild-Type KRAS
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --outdir structures/4ldj/wt

# Example: Create G12D Mutant
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 2.2. Setup and Run MD
**Environment:** `base` | **Directory:** `/content/drive/MyDrive/openmm`

```bash
conda activate base
cd /content/drive/MyDrive/openmm

# 1. Initialize the simulation folder
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1

# 2. Run the pipeline (Example: 1ns)
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 1.0 --traj-interval 10 --equil-time 100 --checkpoint-ps 250
```

> **Modular Control:** You can also run individual steps for more control:
> `colabmda openmm em --name 4ldj_wt`
> `colabmda openmm nvt --name 4ldj_wt --seed 1`
> `colabmda openmm check-equil --name 4ldj_wt`
> `colabmda openmm md --name 4ldj_wt --total-ns 1.0`

> **Note:** The `run` command includes an **Automated Stability Gate**. It automatically analyzes equilibration logs and aborts if the system hasn't stabilized, saving GPU time.

### 2.3. Merge and Center
Combine chunks into a single DCD and wrap solvent.
```bash
# Standard Merge
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap

# Advanced: Save space by keeping every 10th frame
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap --stride 10
```

---

## 📊 3. Analysis & Comparison

### 3.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_wt
```

### 3.2. WT vs Mutant Comparison
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

---

## Project Strategy (WT + Mutants)
Organize work in three phases:
1. **Preparation**: Build WT first in `structures/<pdbid>/wt/`, then generate mutants.
2. **Simulation**: Run WT and mutants in separate folders under `simulations/`.
3. **Analysis**: Store per-system analysis in `analysis/single/`, then generate overlays in `analysis/compare/`.

---

## 📂 Project Structure
```text
/content/drive/MyDrive/openmm/
  structures/
    4ldj/
      wt/          # Wild-type modeled PDBs
      mutants/     # G12D/G12C modeled PDBs
  simulations/
    4ldj_wt/
      r1/          # Replica 1 (em.chk, npt.chk, prod.dcd)
      r2/          # Replica 2
    4ldj_G12D/
      r1/
      r2/
  analysis/
    single/
      4ldj_wt/     # [r1, r2, aggregate] reports
      4ldj_G12D/
    compare/       # Final WT vs Mutant overlays
```

## 📜 Acknowledgements
- **OpenMM & PDBFixer**
- **Modeller**
- **MDAnalysis & MDTraj**
- **NumPy, Matplotlib, Biopython**
- **Google Colab & Miniforge/Conda**

---
**Maintained by**: [Paul Shamrat](https://github.com/paulshamrat)
