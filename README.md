# ColabMDA 🧬
**Publication-grade Molecular Dynamics on Google Colab**

ColabMDA is a modular, high-performance pipeline for protein modeling and MD simulation. It is specifically optimized for Google Colab and Google Drive, with built-in protection against GPU timeouts.

---

## 🚀 Step 1: Installation
These are the **only two commands** you need to get started.

### 1.1. Environment Setup
Install OpenMM, Modeller, and GPU drivers. This takes ~3-5 minutes.
```bash
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest
```

### 1.2. Package Installation
Install the ColabMDA command-line tool.
```bash
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

---

## 🔄 How to Resume After a Timeout
If your Google Colab session expires or the GPU is disconnected:
1. Re-run **Step 1.1 and 1.2** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your checkpoint (`.chk`) files and resume from the last saved state.

---

## 🛠 Step 2: Simulation Workflow

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

## 📊 Step 3: Analysis & Comparison

### 3.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_wt
```

### 3.2. WT vs Mutant Comparison
Generate overlay plots with Mean + Standard Deviation bands.
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

---

## 📝 Troubleshooting & Advanced

<details>
<summary>Fix MODELLER License Key</summary>

If you skipped the key during setup, run this:
```bash
conda activate modeller_env
export KEY_MODELLER='YOUR_KEY'
python3 -c "import os; from colabmda.modeller.utils import patch_modeller; patch_modeller(os.environ.get('KEY_MODELLER'))"
```
</details>

<details>
<summary>Lightweight Installation (OpenMM Only)</summary>

If you don't need Modeller or structural building:
```bash
bash bootstrap_colab_openmm_gpu.sh latest
```
</details>

<details>
<summary>Folder Architecture</summary>

```text
openmm/
  structures/    # Initial PDBs
  simulations/   # MD data (r1, r2)
  analysis/      # Single & Comparative reports
```
</details>

---
**Maintained by**: [Paul Shamrat](https://github.com/paulshamrat)
