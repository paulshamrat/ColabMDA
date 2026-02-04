# ColabMDA

ColabMDA provides Colab-first, terminal-based workflows for molecular modeling and MD: Modeller (CPU) for homology modeling and mutations, and OpenMM (GPU) for protein-in-water simulations. The goal is a clean CLI that wraps stable legacy scripts while staying reliable in short-lived Colab GPU sessions.

Key ideas:
- Colab-first: designed for frequent disconnects and limited storage.
- Terminal-only: no notebooks required for core usage.
- Resume-safe MD: OpenMM runs use checkpoint-based chunks.
- Software-style packaging: required workflow scripts are bundled inside the `colabmda` package.

## Installation (Colab Terminal)

Optional (Colab notebook only): mount Drive and check GPU, then do all installs in the Terminal.

```python
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

> All environment setup and package installation should be performed in the Colab Terminal (not notebook cells).

### 0. Get the Repo

```bash
cd /content
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
```

### 1. Install Miniforge + Mamba

```bash
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
export PATH="$HOME/miniforge3/bin:$PATH"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda install -y -n base -c conda-forge mamba
```

### 2. OpenMM GPU Environment (CUDA)

```bash
mamba create -y -n openmm_gpu -c conda-forge \
  python=3.10 \
  cudatoolkit=11.8 \
  openmm \
  openmmtools \
  pdbfixer \
  mdanalysis \
  mdtraj \
  numpy \
  matplotlib \
  biopython
conda activate openmm_gpu
```

Install ColabMDA into this environment:

```bash
cd /content/ColabMDA
pip install -e .
```

Quick sanity check (optional):

```bash
python - << 'PY'
from openmm import Platform
print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PY
```

### 3. Modeller CPU Environment + License

```bash
conda config --add channels salilab
mamba create -y -n modeller_env python=3.10 modeller biopython
conda activate modeller_env

cd /content/ColabMDA
pip install -e .
```

Set your Modeller license (choose one method):

Option A (recommended for Colab):

```bash
export KEY_MODELLER="YOUR_LICENSE_KEY"
```

Option B (edit config after install):

```bash
LICENSE_KEY="YOUR_LICENSE_KEY"
CONFIG="$HOME/miniforge3/envs/modeller_env/lib/modeller-*/modlib/modeller/config.py"
sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
python - << 'PY'
import modeller
print('Modeller OK, version', modeller.__version__)
PY
```

## Quickstart (CLI)

These commands are safe for Colab and will resume from checkpoints on reconnect.

```bash
conda activate openmm_gpu
cd /content
colabmda openmm prep --pdb-id 4ldj
colabmda openmm run --pdb-id 4ldj --total-ns 1 --traj-interval 100 --checkpoint-ps 100 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj
colabmda openmm status --pdb-id 4ldj
colabmda openmm merge --pdb-id 4ldj
colabmda openmm analysis --pdb-id 4ldj
```

## OpenMM Workflow Overview (CLI)

All outputs go to your **current working directory**. Use `cd` to choose where you want files written. For Colab, use `/content` for speed and stability.

### 1. Prepare the PDB Structure

From RCSB:

```bash
cd /content
colabmda openmm prep --pdb-id 4ldj
```

From a local PDB file:

```bash
colabmda openmm prep --pdb-file /content/4ldj.pdb --name 4ldj --outdir 4ldj
```

Outputs are written to `./<pdbid>/` and include `<pdbid>_cleaned.pdb`.

### 2. Run Chunked MD Simulation

```bash
cd /content
colabmda openmm run --pdb-id 4ldj \
  --total-ns 1 \
  --traj-interval 100 \
  --equil-time 100 \
  --checkpoint-ps 100 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj
```

This runs minimization, equilibration (NVT + NPT), then production in chunks. If the session disconnects, re-run the same command to resume.

### 3. Merge Trajectory Chunks

```bash
colabmda openmm merge --pdb-id 4ldj
```

Outputs `prod_full.dcd` and `prod_full.log` in the same directory.

### 4. Analyze Trajectory

```bash
colabmda openmm analysis --pdb-id 4ldj
```

Outputs go to `analysis_<pdbid>_TIMESTAMP/` unless `--outdir` is set.

### 5. Status Check (Resume Readiness)

```bash
colabmda openmm status --pdb-id 4ldj
```

## Modeller Workflow (CLI)

```bash
colabmda modeller build --pdb-id 4bgq --uniprot-id O76039 --chain A --range 1 303
colabmda modeller mutate --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --mut K76E
```

## Mutant Workflow Walkthrough (WT → Mutant MD)

This shows a full WT + mutant flow using the current CLI. The mutant is built with Modeller, then simulated with OpenMM.

### 1. Prepare WT (OpenMM)

```bash
conda activate openmm_gpu
cd /content
colabmda openmm prep --pdb-id 4ldj
```

### 2. Create Mutant (Modeller)

```bash
conda activate modeller_env
cd /content
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --mut G12C --outdir-mut 4ldj_G12C
```

This produces a mutant PDB at `4ldj_G12C/4ldj_G12C.pdb`.

### 3. Prep Mutant for OpenMM

```bash
conda activate openmm_gpu
cd /content
colabmda openmm prep --pdb-file 4ldj_G12C/4ldj_G12C.pdb --name 4ldj_G12C --outdir 4ldj_G12C
```

### 4. Run Mutant MD (Resume-safe)

```bash
colabmda openmm run --workdir 4ldj_G12C --name 4ldj_G12C --total-ns 1 --traj-interval 100 --checkpoint-ps 100 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj_G12C
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
cd /content
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --list mutations.txt --outdir-mut 4ldj_mutants
```

This creates one mutant PDB per line in `4ldj_mutants/`. For each mutant, run the same OpenMM prep/run/merge/analysis steps as shown above.

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
