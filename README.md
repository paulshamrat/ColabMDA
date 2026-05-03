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
1. Re-run **Steps 1.2 and 1.3** to reinstall the environment.
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

### Modeller License Fix
```bash
conda activate modeller_env
export KEY_MODELLER='YOUR_KEY'
python3 -c "import os; from colabmda.modeller.utils import patch_modeller; patch_modeller(os.environ.get('KEY_MODELLER'))"
```
</details>

---

## 🛠 2. Simulation Workflow

### 2.1. Build Structures (WT and Mutants)
```bash
# Example: Build Wild-Type KRAS
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --outdir structures/4ldj/wt

# Example: Create G12D Mutant
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 2.2. Stage and Run MD
The `run` command handles everything: Minimization, Equilibration (NVT/NPT), Stability Checks, and Production MD.

```bash
# 1. Stage the folder
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1

# 2. Run the pipeline (Example: 1ns)
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 1.0 --traj-interval 10 --equil-time 100
```

### 2.3. Merge and Center
```bash
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap
```

---

## 📊 3. Analysis & Comparison

### 3.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_wt
```

### 3.2. WT vs Mutant Comparison
Generate overlay plots with Mean + Standard Deviation bands.

**Example (Your KRAS Project):**
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

---

## 📂 Project Structure
```text
openmm/
  structures/
    4ldj/
      wt/          # Wild-type modeled PDBs
      mutants/     # G12D/G12C modeled PDBs
  simulations/
    4ldj_wt/       # Wild-type replicas (r1, r2)
    4ldj_G12D/     # Mutant replicas (r1, r2)
  analysis/
    single/        # Individual system results
    compare/       # Aggregate WT vs Mutant overlays
```

## 📜 Acknowledgements
- **OpenMM & PDBFixer**
- **Modeller**
- **MDAnalysis & MDTraj**
- **NumPy, Matplotlib, Biopython**
- **Google Colab & Miniconda/Miniforge**

---
**Maintained by**: [Paul Shamrat](https://github.com/paulshamrat)
