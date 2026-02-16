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

### HPC Quick Install (module-based)

```bash
module purge
module load miniforge3/24.3.0-0
```

Then use the same OpenMM/Modeller env steps above.

## Workflow Overview (WT First, Recommended)

Canonical protocol is:
1. Build/confirm WT with Modeller (`--pdb-id` + `--uniprot-id`)
2. Generate mutants from that WT
3. Run OpenMM for WT and each mutant
4. Compare WT vs mutant analyses

Default project root is `/content/drive/MyDrive/openmm`.

### Folder Layout

```text
/content/drive/MyDrive/openmm/
  4ldj/
    wt/
      modeller/
      openmm/
    mutants/
      4ldj_G12C/
        modeller/
        openmm/
      4ldj_G12D/
        modeller/
        openmm/
```

### 1. Build WT (Modeller, First Step)

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env
cd /content/drive/MyDrive/openmm

colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --outdir 4ldj/wt/modeller
```

Use the produced WT PDB as the input for OpenMM WT prep.

### 2. Create Mutants from WT

```bash
# Single mutant
colabmda modeller mutate --pdb-in 4ldj/wt/modeller/<wt_model>.pdb --chain A --mut G12C --outdir-mut 4ldj/mutants/4ldj_G12C/modeller

# Batch mutants
colabmda modeller mutate --pdb-in 4ldj/wt/modeller/<wt_model>.pdb --chain A --list mutations.txt --outdir-mut 4ldj/mutants
```

### 3. Run OpenMM for WT

```bash
conda activate base

# Prep from WT model file
colabmda openmm prep --pdb-file 4ldj/wt/modeller/<wt_model>.pdb --name 4ldj_wt --outdir 4ldj/wt/openmm/run

# Run / merge / analysis
colabmda openmm run --workdir /content/drive/MyDrive/openmm/4ldj/wt/openmm/run --name 4ldj_wt --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --pdb-dir /content/drive/MyDrive/openmm/4ldj/wt/openmm/run --stride 10
colabmda openmm analysis --pdb-dir /content/drive/MyDrive/openmm/4ldj/wt/openmm/run --interval 10 --outdir /content/drive/MyDrive/openmm/4ldj/wt/openmm/analysis
```

### 4. Run OpenMM for Each Mutant

```bash
# Example: G12C
colabmda openmm prep --pdb-file 4ldj/mutants/4ldj_G12C/modeller/4ldj_G12C.pdb --name 4ldj_G12C --outdir 4ldj/mutants/4ldj_G12C/openmm/run
colabmda openmm run --workdir /content/drive/MyDrive/openmm/4ldj/mutants/4ldj_G12C/openmm/run --name 4ldj_G12C --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --pdb-dir /content/drive/MyDrive/openmm/4ldj/mutants/4ldj_G12C/openmm/run --stride 10
colabmda openmm analysis --pdb-dir /content/drive/MyDrive/openmm/4ldj/mutants/4ldj_G12C/openmm/run --interval 10 --outdir /content/drive/MyDrive/openmm/4ldj/mutants/4ldj_G12C/openmm/analysis
```

### 5. Quick OpenMM-Only Path (No Modeller)

If you only want direct protein-in-water MD from PDB ID:

```bash
colabmda openmm prep --pdb-id 4ldj
colabmda openmm run --pdb-id 4ldj --total-ns 5 --traj-interval 1 --equil-time 100 --checkpoint-ps 1000
colabmda openmm merge --pdb-id 4ldj --stride 10
colabmda openmm analysis --pdb-id 4ldj --interval 10
```

## Project Strategy (WT + Mutants)

Organize work in three phases:

1. Preparation (WT and mutants)
Place the original PDB in a dedicated folder, then create `wt/` and `mutants/` subfolders. Use Modeller to confirm or revert mutations as needed using UniProt as the reference sequence.
2. Simulation (WT and mutants)
Run OpenMM for `wt/` first, then each mutant in its own folder.
3. Analysis (WT + mutants)
Generate RMSD/RMSF/Rg for each system, then plot WT and all mutants together for comparison.

## Where Files Go

- `openmm prep --pdb-id 4ldj` writes to `/content/drive/MyDrive/openmm/4ldj/prep/`
- `openmm prep --pdb-file ... --name 4ldj_g12c` writes to `/content/drive/MyDrive/openmm/4ldj_g12c/prep/`
- `openmm run --pdb-id 4ldj` writes to `/content/drive/MyDrive/openmm/4ldj/run/`
- `openmm merge --pdb-id 4ldj` writes merged files in `/content/drive/MyDrive/openmm/4ldj/run/`
- `openmm analysis --pdb-id 4ldj` writes plots under `/content/drive/MyDrive/openmm/4ldj/analysis/`

## Colab Limitations and Best Practices

- Default workflow writes directly to Google Drive for persistence.
- For faster local SSD runs, you can still use explicit `--workdir /content/work/...` plus `--sync-dir`.
- Avoid saving trajectories too frequently. Recommended: 1 frame per 100 ps.
- Do not write large trajectories directly to Google Drive.
- Expect GPU disconnects after a few hours; resume-safe runs are mandatory.
- On a new Colab session (same day or next day), re-run the installation steps and then run the exact same `colabmda openmm run` command to resume from checkpoints.

## HPC (Slurm) Workflow

This section mirrors the CLI workflow for HPC systems using modules and Slurm. Adjust module names, GPU types, and time limits to match your cluster.

### 1. Prepare Environments (login node)

```bash
module purge
module load miniforge3/24.3.0-0
cd /path/to/ColabMDA
bash envs/install_openmm_env.sh
```

### 2. Recommended Runtime Settings

```bash
export OPENMM_DEFAULT_PLATFORM=CUDA
export OPENMM_CUDA_DEFAULT_PRECISION=mixed
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
```

Sanity check:

```bash
python - << 'PY'
from openmm import Platform
print("Platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PY
```

### 3. Batch Script Example (Slurm)

Save as `run_on_hpc.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=4bgq_10ns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=a100:1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00

module purge
module load miniforge3/24.3.0-0
source activate pw310

PDBID="4bgq"
TOTAL_NS=5.0
INTERVAL_PS=1.0
EQUIL_PS=100.0
CHUNK_PS=1000.0

export OPENMM_DEFAULT_PLATFORM=CUDA
export OPENMM_CUDA_DEFAULT_PRECISION=mixed
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

echo "Step 1: Cleaning PDB..."
colabmda openmm prep --pdb-id "$PDBID"

echo "Step 2: Running Production..."
colabmda openmm run --pdb-id "$PDBID" \
  --total-ns "$TOTAL_NS" \
  --traj-interval "$INTERVAL_PS" \
  --equil-time "$EQUIL_PS" \
  --checkpoint-ps "$CHUNK_PS"

echo "Step 3: Merging Trajectories..."
colabmda openmm merge --pdb-id "$PDBID"

echo "Step 4: Analysis..."
colabmda openmm analysis --pdb-id "$PDBID"

echo "Done."
```

Submit with:

```bash
sbatch run_on_hpc.sh
```

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
