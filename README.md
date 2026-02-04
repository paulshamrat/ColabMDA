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

### 2. Installation on Terminal

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

### 3. Install ColabMDA

```bash
cd /content/drive/MyDrive/openmm_runs
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
pip install -e .
```

### 4. Modeller CPU Environment + License

```bash
conda config --add channels salilab
mamba create -y -n modeller_env python=3.10 modeller biopython
conda activate modeller_env

cd /content/drive/MyDrive/openmm_runs/ColabMDA
pip install -e .
```

Set your Modeller license (required to run Modeller):

```bash
export KEY_MODELLER="MODELIRANJE"
```

### HPC Quick Install (module-based)

```bash
module purge
module load miniforge3/24.3.0-0
```

Then use the same OpenMM/Modeller env steps above.

## Workflow Overview (Updated CLI)

All commands run from your working directory. For Colab, use Drive so outputs persist:

```bash
cd /content/drive/MyDrive/openmm_runs
```

### 1. Prepare the PDB Structure

Download and clean a PDB structure, removing heterogens (except water), building missing residues/atoms, and adding hydrogens at pH 7.0.

```bash
colabmda openmm prep --pdb-id 4ldj
```

Output: creates `4ldj/` containing `4ldj.pdb` (raw) and `4ldj_cleaned.pdb` (processed).

### 2. Run Chunked MD Simulation (5 ns, 1 ps frames, 1000 ps chunks)

```bash
colabmda openmm run --pdb-id 4ldj \
  --total-ns 5 \
  --traj-interval 1 \
  --equil-time 100 \
  --checkpoint-ps 1000
```

Output: multiple chunk files (`prod_*`), logs, checkpoints, and system files in `4ldj/`.

### 3. Merge Trajectory Chunks

```bash
colabmda openmm merge --pdb-id 4ldj
```

Output: merged trajectory (`prod_full.dcd`) and log (`prod_full.log`).

### 4. Analyze Trajectory

```bash
colabmda openmm analysis --pdb-id 4ldj
```

Output: analysis results in `analysis_<pdbid>_TIMESTAMP/`.

## Modeller Workflow (CLI)

```bash
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169
```

Before making mutants, confirm whether the PDB is WT or already mutated by comparing against the UniProt reference sequence (P01116). If the PDB is already a mutant (e.g., G12C), flag it and **revert to WT first**, then generate mutants from the WT.

Recommended folder layout under Drive:

```
/content/drive/MyDrive/openmm_runs/4ldj/
  original/
  wt/
  mutants/
```

Example (revert mutant to WT, then make mutants):

```bash
# Revert G12C -> WT (C12G)
colabmda modeller mutate --pdb-in 4ldj/original/4ldj_cleaned.pdb --chain A --mut C12G --outdir-mut 4ldj/wt

# Now make mutants from the WT
colabmda modeller mutate --pdb-in 4ldj/wt/4ldj_wt.pdb --chain A --mut G12C --outdir-mut 4ldj/mutants/4ldj_G12C
```

Example mutant list (one per line):

```text
G12C
I36M
G60R
T58I
```

## Mutant Workflow Walkthrough (WT → Mutant MD)

This shows a full WT + mutant flow using the current CLI. The mutant is built with Modeller, then simulated with OpenMM.

### 1. Prepare WT (OpenMM)

```bash
conda activate openmm_gpu
cd /content/drive/MyDrive/openmm_runs   # or: cd /content
colabmda openmm prep --pdb-id 4ldj
```

### 2. Create Mutant (Modeller)

```bash
conda activate modeller_env
cd /content/drive/MyDrive/openmm_runs   # or: cd /content
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --mut G12C --outdir-mut 4ldj_G12C
```

This produces a mutant PDB at `4ldj_G12C/4ldj_G12C.pdb`.

### 3. Prep Mutant for OpenMM

```bash
conda activate openmm_gpu
cd /content/drive/MyDrive/openmm_runs   # or: cd /content
colabmda openmm prep --pdb-file 4ldj_G12C/4ldj_G12C.pdb --name 4ldj_G12C --outdir 4ldj_G12C
```

### 4. Run Mutant MD (Resume-safe)

```bash
colabmda openmm run --workdir 4ldj_G12C --name 4ldj_G12C --total-ns 5 --traj-interval 1 --checkpoint-ps 1000
colabmda openmm status --pdb-id 4ldj_G12C
colabmda openmm merge --pdb-dir 4ldj_G12C
colabmda openmm analysis --pdb-dir 4ldj_G12C
```

## Batch Mutants (List File)

If you have many mutations, put one per line in a text file:

```text
G12C
G12D
Q61L
```

Then run Modeller in batch mode:

```bash
conda activate modeller_env
cd /content/drive/MyDrive/openmm_runs   # or: cd /content
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --list mutations.txt --outdir-mut 4ldj_mutants
```

This creates one mutant PDB per line in `4ldj_mutants/`. For each mutant, run the same OpenMM prep/run/merge/analysis steps as shown above.

## Project Strategy (WT + Mutants)

Organize work in three phases:

1. Preparation (WT and mutants)
Place the original PDB in a dedicated folder, then create `wt/` and `mutants/` subfolders. Use Modeller to confirm or revert mutations as needed using UniProt as the reference sequence.
2. Simulation (WT and mutants)
Run OpenMM for `wt/` first, then each mutant in its own folder.
3. Analysis (WT + mutants)
Generate RMSD/RMSF/Rg for each system, then plot WT and all mutants together for comparison.

## Where Files Go

- `openmm prep --pdb-id 4ldj` writes to `./4ldj/`
- `openmm prep --pdb-file ... --outdir 4ldj_g12c` writes to `./4ldj_g12c/`
- `openmm run` writes trajectories, checkpoints, logs inside the work directory
- `openmm merge` writes merged files in the same work directory
- `openmm analysis` writes plots to `analysis_<pdbid>_TIMESTAMP/` unless `--outdir` is set

## Colab Limitations and Best Practices

- Run MD in `/content/` (local SSD) for speed and stability.
- Sync only completed chunks to Google Drive.
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
