# ColabMDA


**ColabMDA** is a specialized tool that lets you run high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your work. The biggest problem with Colab is that it disconnects, often destroying hours of simulation data. **ColabMDA fixes this** with a "resume-safe" system that automatically saves your progress to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From modeling protein mutations to generating publication-ready analysis, ColabMDA handles the complex setup for you, making it the easiest way to get high-quality MD results using free cloud GPUs.

📖 **Full Documentation:** Visit our official manual at [colabmda.readthedocs.io](https://colabmda.readthedocs.io/)

| Category | Details |
| :--- | :--- |
| **Release** | [![GitHub tag](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) |
| **Availability** | [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) |
| **Documentation** | [![Documentation Status](https://readthedocs.org/projects/colabmda/badge/?version=latest)](https://colabmda.readthedocs.io/en/latest/?badge=latest) |
| **Workflows** | [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) |
| **Issues** | [![GitHub issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) |
| **License** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) |
| **Style / Lint** | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |
| **Dependencies** | `OpenMM`, `Modeller`, `MDAnalysis`, `MDTraj` |
| **Platform** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/05-colabmd-simulation-2024.ipynb) `Linux` `HPC (SLURM)` |

---

## 1. Installation

> 💡 **Terminal Access:** All bash commands should be run in the **Colab Terminal** (Open via the **⋮** menu -> **Terminal**).

### 1.1. Setup Colab Runtime & Drive
Before starting, ensure your environment is ready:
1. **Enable GPU:** Go to `Runtime` -> `Change runtime type` and select **T4 GPU**.
2. **Verify GPU:** Run `!nvidia-smi` in a cell to confirm GPU access.
3. **Mount Drive:** Click the **Folder icon** 📂 in the left sidebar, then click the **Drive icon** (Mount Drive), or run the code block below:

```python
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

### 1.2. Environment & Package Installation (Required)
Run the following in the **Colab Terminal** (⋮ → Terminal). 
*Estimated time: ~3–5 minutes.*

```bash
# 1. Install the core scientific environment
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest

# 2. Install ColabMDA package
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

> [!IMPORTANT]
> **Modeller License Prompt:** During Step 1, the script will **pause** and ask you to `Enter your Modeller License Key`. You must paste your key and press **Enter** to proceed. The installation will not complete without it.
>
> 🔑 **Get a Free License:** If you don't have one, register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Academic licenses are free and sent instantly via email).

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Steps 1.2** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your `.chk` files and resume from where it left off.

---

<details>
<summary>Manual / Detailed Installation (Advanced)</summary>

### A. Manual Terminal Installation (Step-by-Step)
In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
# Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

# Step 2: Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Step 3: Install Mamba into the base environment
conda install -y -n base -c conda-forge mamba

# Step 4: Install CUDA-enabled OpenMM and OpenMMTools
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools

# Step 5: Install PDBFixer (conda, fallback to pip)
conda install -y -c conda-forge pdbfixer || pip install pdbfixer

# Step 6: Install MDAnalysis, MDTraj, NumPy, Matplotlib, and Biopython
mamba install -y -c conda-forge mdanalysis mdtraj numpy matplotlib biopython

# Step 7: Verify installations
python3 - << 'EOF'
from openmm import Platform; print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import MDAnalysis, mdtraj, Bio; print("MDAnalysis:", MDAnalysis.__version__, "MDTraj:", mdtraj.__version__, "Biopython:", Bio.__version__)
EOF
```

### B. Alternative: Script-based Installation
```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/install_colabmda_release.sh -o install_colabmda_release.sh
bash install_colabmda_release.sh latest /content/colabmda
```

### C. Modeller CPU Environment Setup
```bash
cd /content/drive/MyDrive/openmm/ColabMDA
bash envs/install_modeller_env.sh
```
</details>

---

## 2. Beyond Colab: Local & HPC Usage
ColabMDA is not limited to the cloud. It is a full-featured MD pipeline that works on any Linux system with an NVIDIA GPU.

### Local Workstation Setup
Use the provided `environment.yml` to create a production-ready environment in one command:
```bash
mamba env create -f environment.yml
conda activate colabmda
```

### HPC Usage (SLURM)
You can easily incorporate ColabMDA into SLURM batch scripts. Since it processes trajectories in chunks, it is highly efficient for long-running jobs on cluster partitions with time limits.

---

## 3. Simulation Workflow

### Pipeline at a Glance:
1. **Model:** Build structures (WT & Mutants) in `structures/`
2. **Stage:** Initialize the simulation folder in `simulations/`
3. **Run:** Execute the MD simulation (Resume-safe)
4. **Merge:** Combine trajectory chunks into a final file
5. **Analyze:** Generate RMSD, Rg, and RMSF plots

### 3.1. Build Structures (WT and Mutants)
**Environment:** `modeller_env`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env
cd /path/to/your/project

# Example: Build Wild-Type KRAS
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --outdir structures/4ldj/wt

# Example: Create G12D Mutant
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 3.2. Setup and Run MD
**Environment:** `base`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate base
cd /path/to/your/project

# 1. Initialize the simulation folder
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1

# 2. Run the pipeline (Example: 5ns)
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 5.0 --traj-interval 1 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Storage Tip:** For a typical system (e.g., KRAS in water, ~30,000 atoms), a 100ns run at high resolution (1ps) can produce over **36GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~3.6GB. Always calculate your storage needs based on your specific system size before starting long runs.

> **Note:** The `run` command includes an **Automated Stability Gate**. It automatically analyzes equilibration logs and aborts if the system hasn't stabilized, saving GPU time.

### 3.3. Merge and Center
Combine chunks into a single DCD and wrap solvent.

```bash
# Standard Merge (Center + Wrap)
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap
```

> 💡 **Pro-Tip for Long Runs:**
> Merging processes trajectories frame-by-frame, so it won't crash your RAM. You can merge without striding (`--stride 1`) for full resolution, or use `--stride 10` to create a lightweight file for local viewing.

---

## 4. Analysis & Comparison

### 4.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_wt
```

> ⚠️ **Analysis Tip:** If your plots show the wrong time scale (e.g., 10ns instead of 100ns), provide the frame interval manually. For example, if you ran with `--traj-interval 10`:
> `colabmda openmm analysis --pdb-id 4ldj_wt --interval 10`

### 4.2. WT vs Mutant Comparison
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

---

## Project Strategy
Organize work in three phases:
1. **Preparation**: Build WT first in `structures/<pdbid>/wt/`, then generate mutants.
2. **Simulation**: Run WT and mutants in separate folders under `simulations/`.
3. **Analysis**: Store per-system analysis in `analysis/single/`, then generate overlays in `analysis/compare/`.

---

## Project Structure
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

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

## Changelog

### v0.1.0 (Initial Beta)
*   **Modular Pipeline:** New modular CLI for EM, NVT, NPT, and Production MD.
*   **Resume-Safe Engine:** Integrated checkpointing logic for fail-safe simulations on Google Colab.
*   **Modeling:** Automated Wild-Type building and mutation support via Modeller.
*   **Analysis:** Robust trajectory merging and comparative RMSD/Rg/RMSF analysis tools.
*   **Professional Standards:** Added CI/CD workflows, Black formatting, and Ruff linting.

---

## Citation

This repository was inspired by the methodologies established in the research published below. Originally developed as a simple GROMACS-on-Colab workflow, ColabMDA has since evolved into a specialized OpenMM-centered pipeline. If you use this tool, please consider citing the underlying study:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
