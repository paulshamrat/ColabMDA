# ColabMDA

ColabMDA provides Colab-first, terminal-based workflows for molecular modeling and MD: Modeller (CPU) for homology modeling and mutations, and OpenMM (GPU) for protein-in-water simulations. The goal is a clean CLI that wraps stable legacy scripts while staying reliable in short-lived Colab GPU sessions.

Key ideas:
- Colab-first: designed for frequent disconnects and limited storage.
- Terminal-only: no notebooks required for core usage.
- Resume-safe MD: OpenMM runs use checkpoint-based chunks.
- Software-style packaging: required workflow scripts are bundled inside the `colabmda` package.

## 1. Installation

<details>
<summary><h3>1.1. Notebook Setup (Colab)</h3></summary>

Open a new Colab notebook in Google Drive. In the first cell, mount Drive and confirm GPU access:

```python
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

> All environment setup and package installation (Conda, Mamba, OpenMM, analysis libraries) should be performed in the Colab Terminal, not notebook cells.
</details>

<details open>
<summary><h3>1.2. Environment Setup (Required)</h3></summary>

In the Colab Terminal, run this to install **OpenMM + Modeller** in one flow (recommended):

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest
```

```bash
bash bootstrap_colab_openmm_gpu.sh latest
```
</details>

### 1.3. Package Installation (Required)
```bash
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Steps 1.2 and 1.3** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your `.chk` files and resume from where it left off.

<details>
<summary><h3>1.3. Installation on Terminal (Manual Alternative)</h3></summary>

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
</details>

<details>
<summary><h3>1.4. Install ColabMDA (No Full Clone)</h3></summary>

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
</details>

<details>
<summary><h3>1.5. Modeller CPU Environment + License</h3></summary>

```bash
cd /content/drive/MyDrive/openmm/ColabMDA
bash envs/install_modeller_env.sh
```

The script persists the key for the conda environment/future shells and patches MODELLER `config.py` automatically.

> **Note:** If the interactive prompt is skipped in Colab (you see a warning), you can set the key manually after installation:
> ```bash
> conda activate modeller_env
> export KEY_MODELLER='YOUR_KEY'
> python3 -c "import os; from colabmda.modeller.utils import patch_modeller; patch_modeller(os.environ.get('KEY_MODELLER'))"
> ```
</details>


## 2. Simulation Workflow

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
      r1/
      r2/
    4ldj_G12D/
      r1/
      r2/
  analysis/
    single/
      4ldj_wt/
        r1/
        r2/
        aggregate/
      4ldj_G12D/
        r1/
        r2/
        aggregate/
    compare/
      wt_vs_mutants/
```

### 2.1. Build WT and Mutants in `structures/`

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

### 2.2. Stage One Structure into `simulations/<name>`

```bash
conda activate base
cd /content/drive/MyDrive/openmm

# Stage WT (Replica 1)
colabmda openmm stage --pdb-file structures/4ldj/wt/<wt_model>.pdb --name 4ldj_wt --replica r1
```
**Example:**
```bash
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1
```

# Stage one mutant (Replica 1)
colabmda openmm stage --pdb-file structures/4ldj/mutants/4ldj_G12C/4ldj_G12C.pdb --name 4ldj_G12C --replica r1

This creates:
- `/content/drive/MyDrive/openmm/simulations/4ldj_wt/r1`
- `/content/drive/MyDrive/openmm/simulations/4ldj_G12C/r1`

### 2.3. Run from Inside Simulation Folder

```bash
cd /content/drive/MyDrive/openmm/simulations/4ldj_wt/r1

# One-click Modular Run (EM -> NVT -> NPT -> Stability Check -> MD)
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000 --seed 1

# Merge with automatic PBC centering and solvent wrapping
colabmda openmm merge --stride 10 --center --wrap

## 3. Analysis & Comparison

### 3.1. Standard Analysis
colabmda openmm analysis --pdb-id 4ldj_wt
```

> **Modular Control:** You can also run individual steps for more control:
> `colabmda openmm em --name 4ldj_wt`
> `colabmda openmm nvt --name 4ldj_wt --seed 1`
> `colabmda openmm check-equil --name 4ldj_wt`
> `colabmda openmm md --name 4ldj_wt --total-ns 5`

> **Note:** The `run` command now includes an **Automated Stability Gate**. It will automatically analyze your equilibration logs and abort if the system hasn't stabilized, saving you valuable GPU time.

For mutant:

```bash
cd /content/drive/MyDrive/openmm
colabmda openmm run --name 4ldj_G12D --replica r1 --total-ns 1.0 --traj-interval 10 --equil-time 100
colabmda openmm merge --pdb-dir simulations/4ldj_G12D/r1 --center --wrap
colabmda openmm analysis --pdb-id 4ldj_G12D
```

### 3.2. Compare WT vs Mutants (Aggregate Overlay)

Once you have analyzed both WT and Mutants, generate the final publication comparison:

```bash
colabmda openmm compare \
  --series "LABEL=DIR1,DIR2" \
  --series "LABEL2=DIR3,DIR4" \
  --outdir analysis/compare/project_name
```

**Example (Your KRAS Project):**
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

## 4. Advanced Options

### 4.1. Quick OpenMM-Only Path (No Modeller)

```bash
colabmda openmm prep --pdb-id 4ldj
colabmda openmm run --pdb-id 4ldj --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --pdb-id 4ldj --stride 10
colabmda openmm analysis --pdb-id 4ldj
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
- Staged simulation folders: `/content/drive/MyDrive/openmm/simulations/<name>/<replica>/`
- Single-system analysis: `/content/drive/MyDrive/openmm/analysis/single/<name>/[r1, r2, aggregate]/`
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
