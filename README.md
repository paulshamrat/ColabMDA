# ColabMDA

ColabMDA provides Colab-first, terminal-based workflows for molecular modeling and MD: Modeller (CPU) for homology modeling and mutations, and OpenMM (GPU) for protein-in-water simulations. The goal is a clean CLI that wraps stable legacy scripts while staying reliable in short-lived Colab GPU sessions.

Key ideas:
- Colab-first: designed for frequent disconnects and limited storage.
- Terminal-only: no notebooks required for core usage.
- Resume-safe MD: OpenMM runs use checkpoint-based chunks.
- Software-style packaging: required workflow scripts are bundled inside the `colabmda` package.

## Installation

### 1. Notebook Setup (Colab)

Open a new Colab notebook in Google Drive. In the first cell, mount Drive and confirm GPU access:

```python
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

> All environment setup and package installation (Conda, Mamba, OpenMM, analysis libraries) should be performed in the Colab Terminal, not notebook cells.

### 2. One-Command Colab Bootstrap (Recommended)

In the Colab Terminal, run:

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
bash bootstrap_colab_openmm_gpu.sh latest
```

Install OpenMM + Modeller in one flow (with MODELLER license prompt):

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest
```

This will:
- install Miniforge (if missing)
- install OpenMM + analysis stack in conda `base` (legacy-compatible path)
- optionally install MODELLER in `modeller_env` when `WITH_MODELLER=1`
- install `colabmda` from GitHub Release wheel (no full repo clone)
- validate GPU/OpenMM platforms and CLI
- create `/content/work` and `/content/drive/MyDrive/openmm`

### 3. Installation on Terminal (Manual Alternative)

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
## 00 installation
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

### 4. Install ColabMDA (No Full Clone)

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/install_colabmda_release.sh -o install_colabmda_release.sh
bash install_colabmda_release.sh latest /content/colabmda
```

This installs from the latest GitHub Release wheel and avoids cloning the full repository.

Developer-only option (if you need source editing):

```bash
cd /content/drive/MyDrive/openmm
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
pip install -e .
```

### 5. Modeller CPU Environment + License

```bash
cd /content/drive/MyDrive/openmm/ColabMDA
bash envs/install_modeller_env.sh
```

During installation, the script prompts for MODELLER license key (`KEY_MODELLER`).
Example key:

```bash
MODELIRANJE
```

The script persists the key for the conda environment and future shells.


## Workflow Overview (WT First, Recommended)

Canonical protocol:
1. Build WT and mutant structures in one place (`structures/`)
2. Stage one selected structure into simulation folder (`simulations/<name>`)
3. `cd` into that folder and run short commands (no long paths)
4. Save per-system analysis in `analysis/single/` and compare in `analysis/compare/`

Default project root: `/content/drive/MyDrive/openmm`

### Folder Layout

```text
/content/drive/MyDrive/openmm/
  structures/
    4ldj/
      wt/
      mutants/
  simulations/
    4ldj_wt/
    4ldj_G12C/
    4ldj_G12D/
  analysis/
    single/
    compare/
```

### 1. Build WT and Mutants in `structures/`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env
cd /content/drive/MyDrive/openmm

# WT
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --outdir structures/4ldj/wt

# Mutants from WT
colabmda modeller mutate --pdb-in structures/4ldj/wt/<wt_model>.pdb --chain A --mut G12C --outdir-mut structures/4ldj/mutants/4ldj_G12C
colabmda modeller mutate --pdb-in structures/4ldj/wt/<wt_model>.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 2. Stage One Structure into `simulations/<name>`

```bash
conda activate base
cd /content/drive/MyDrive/openmm

# Stage WT
colabmda openmm stage --pdb-file structures/4ldj/wt/<wt_model>.pdb --name 4ldj_wt

# Stage one mutant
colabmda openmm stage --pdb-file structures/4ldj/mutants/4ldj_G12C/4ldj_G12C.pdb --name 4ldj_G12C
```

This creates:
- `/content/drive/MyDrive/openmm/simulations/4ldj_wt`
- `/content/drive/MyDrive/openmm/simulations/4ldj_G12C`

### 3. Run from Inside Simulation Folder (Short Commands)

```bash
cd /content/drive/MyDrive/openmm/simulations/4ldj_wt
colabmda openmm run --name 4ldj_wt --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --stride 10
colabmda openmm analysis --interval 10 --outdir /content/drive/MyDrive/openmm/analysis/single/4ldj_wt
```

For mutant:

```bash
cd /content/drive/MyDrive/openmm/simulations/4ldj_G12C
colabmda openmm run --name 4ldj_G12C --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --stride 10
colabmda openmm analysis --interval 10 --outdir /content/drive/MyDrive/openmm/analysis/single/4ldj_G12C
```

### 4. Compare WT vs Mutants in `analysis/compare`

```bash
python openmm/openmm_proteinwater/openmm_compare_plots.py \
  --series WT=/content/drive/MyDrive/openmm/analysis/single/4ldj_wt \
  --series G12C=/content/drive/MyDrive/openmm/analysis/single/4ldj_G12C \
  --series G12D=/content/drive/MyDrive/openmm/analysis/single/4ldj_G12D \
  --outdir /content/drive/MyDrive/openmm/analysis/compare/4ldj_wt_vs_mutants
```

### 5. Quick OpenMM-Only Path (No Modeller)

```bash
colabmda openmm prep --pdb-id 4ldj
colabmda openmm run --pdb-id 4ldj --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --pdb-id 4ldj --stride 10
colabmda openmm analysis --pdb-id 4ldj --interval 10
```

## Project Strategy (WT + Mutants)

Organize work in three phases:

1. Preparation (WT and mutants)
Build WT first in `structures/<pdbid>/wt/`, then generate mutants in `structures/<pdbid>/mutants/`.
2. Simulation (WT and mutants)
Run WT and mutants in separate folders under `simulations/` (for example: `simulations/4ldj_wt`, `simulations/4ldj_G12C`).
3. Analysis (WT + mutants)
Store per-system analysis in `analysis/single/`, then generate overlays in `analysis/compare/`.

## Where Files Go

- WT/mutant structures: `/content/drive/MyDrive/openmm/structures/...`
- Staged simulation folders: `/content/drive/MyDrive/openmm/simulations/<name>/`
- Single-system analysis: `/content/drive/MyDrive/openmm/analysis/single/<name>/`
- WT vs mutant overlays: `/content/drive/MyDrive/openmm/analysis/compare/<project>/`

## Colab Limitations and Best Practices

- Default workflow writes directly to Google Drive for persistence.
- For faster local SSD runs, you can still use explicit `--workdir /content/work/...` and optionally `--sync-dir`.
- Avoid saving trajectories too frequently. Recommended: 1 frame per 100 ps.
- Do not write large trajectories directly to Google Drive.
- Expect GPU disconnects after a few hours; resume-safe runs are mandatory.
- On a new Colab session (same day or next day), re-run the installation steps and then run the exact same `colabmda openmm run` command to resume from checkpoints.

## Notes

- This repo preserves legacy research scripts, but the CLI runs bundled copies inside the package for reliability.
- Keep dependencies minimal in `pyproject.toml`; OpenMM and Modeller are installed via Conda environments.

## Acknowledgements

- OpenMM
- PDBFixer
- MDAnalysis
- MDTraj
- NumPy
- Matplotlib
- Biopython
- Google Colab
- Miniforge, Conda, Mamba

For questions or issues, please contact the repository maintainer.
