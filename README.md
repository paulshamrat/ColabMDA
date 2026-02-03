# ColabMDA

ColabMDA provides Colab-first, terminal-based workflows for molecular modeling and MD: Modeller (CPU) for homology modeling and mutations, and OpenMM (GPU) for protein-in-water simulations. The goal is a clean CLI that wraps stable legacy scripts while staying reliable in short-lived Colab GPU sessions.

Key ideas:
- Colab-first: designed for frequent disconnects and limited storage.
- Terminal-only: no notebooks required.
- Resume-safe MD: OpenMM runs use checkpoint-based chunks.
- Software-style packaging: required workflow scripts are bundled inside the `colabmda` package.

## Installation (Colab Terminal)

All setup should be done in the Colab Terminal (not notebook cells). These steps also work on any Linux machine.
Start in `/content` so all outputs land on Colab's local SSD:

```bash
cd /content
```

### 1. Install Miniforge

```bash
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
bash /tmp/miniforge.sh -b -p "$HOME/miniforge"
export PATH="$HOME/miniforge/bin:$PATH"
source "$HOME/miniforge/etc/profile.d/conda.sh"
```

### 2. OpenMM GPU Environment (CUDA)

Create an OpenMM environment with CUDA enabled:

```bash
conda create -y -n openmm_gpu -c conda-forge \
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

# Install ColabMDA into this env
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

Modeller is **not** on conda-forge. Install it from the official `salilab` channel:

```bash
conda config --add channels salilab
conda create -y -n modeller_env python=3.10 modeller biopython
conda activate modeller_env

# Install ColabMDA into this env (so the CLI is available)
cd /content/ColabMDA
pip install -e .
```

Get a free academic license at the Modeller site, then set it (choose one method):

Option A: set `KEY_MODELLER` before install (recommended for Colab):

```bash
export KEY_MODELLER="YOUR_LICENSE_KEY"
```

Option B: edit the config after install:

```bash
LICENSE_KEY="YOUR_LICENSE_KEY"
CONFIG="$HOME/miniforge/envs/modeller_env/lib/modeller-*/modlib/modeller/config.py"
sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
python - << 'PY'
import modeller
print('Modeller OK, version', modeller.__version__)
PY
```

## Install ColabMDA (Editable)

Clone the repo, then install from the repository root:

```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
pip install -e .
```

Note: the dot at the end is required. `pip install -e` without `.` will fail.

## Workflow (CD-Based, Step-by-Step)

All outputs go to your **current working directory**. Use `cd` to choose where you want files written.

### 1) Install (once per environment)

OpenMM environment:

```bash
conda activate openmm_gpu
cd /content/ColabMDA
pip install -e .
```

Modeller environment:

```bash
conda activate modeller_env
cd /content/ColabMDA
pip install -e .
```

### 2) Prepare Structures (WT first, then mutants)

Work in `/content` for speed and stability:

```bash
cd /content
```

WT prep:

```bash
colabmda openmm prep --pdb-id 4ldj
```

Mutant prep (run in `modeller_env`):

```bash
conda activate modeller_env
cd /content
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --mut G12C --outdir-mut 4ldj_G12C
```

### 3) Simulate WT, then Mutants (OpenMM GPU)

WT simulation:

```bash
conda activate openmm_gpu
cd /content
colabmda openmm run --pdb-id 4ldj --total-ns 1 --traj-interval 100 --checkpoint-ps 100 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj
colabmda openmm status --pdb-id 4ldj
colabmda openmm merge --pdb-id 4ldj
colabmda openmm analysis --pdb-id 4ldj
```

Mutant simulation:

```bash
conda activate openmm_gpu
cd /content
colabmda openmm prep --pdb-file 4ldj_G12C/4ldj_G12C.pdb --name 4ldj_G12C --outdir 4ldj_G12C
colabmda openmm run --workdir 4ldj_G12C --name 4ldj_G12C --total-ns 1 --traj-interval 100 --checkpoint-ps 100 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj_G12C
colabmda openmm status --pdb-id 4ldj_G12C
colabmda openmm merge --pdb-dir 4ldj_G12C
colabmda openmm analysis --pdb-dir 4ldj_G12C
```

### Where Each Command Writes Files

- `openmm prep --pdb-id 4ldj` → `./4ldj/`
- `openmm prep --pdb-file ... --outdir 4ldj_g12c` → `./4ldj_g12c/`
- `openmm run` → all trajectories, checkpoints, logs inside the work directory
- `openmm merge` → writes merged files in the same work directory
- `openmm analysis` → writes plots to `analysis_<pdbid>_TIMESTAMP/` unless `--outdir` is set

To keep data across Colab disconnects, always use `--sync-dir` with `openmm run`. This copies essential outputs to Drive after each chunk. If a runtime disconnects, re-run the same command and it will resume from the checkpoint in the synced folder.

### Modeller (CPU) Workflow

```bash
colabmda modeller build --pdb-id 4bgq --uniprot-id O76039 --chain A --range 1 303
colabmda modeller mutate --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --mut K76E
```

## Colab Limitations and Best Practices

- Run MD in `/content/` (local SSD) for speed and stability.
- Sync only completed chunks to Google Drive.
- Avoid saving trajectories too frequently. Recommended: 1 frame per 100 ps.
- Do not write large trajectories directly to Google Drive.
- Expect GPU disconnects after a few hours; resume-safe runs are mandatory.

## Notes

- This repo preserves legacy research scripts, but the CLI runs **bundled** copies inside the package for reliability.
- Keep dependencies minimal in `pyproject.toml`; OpenMM and Modeller are installed via Conda environments.

If anything in the setup is unclear, open an issue with your exact terminal output.
