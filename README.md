# ColabMDA

ColabMDA provides Colab-first, terminal-based workflows for molecular modeling and MD: Modeller (CPU) for homology modeling and mutations, and OpenMM (GPU) for protein-in-water simulations. The goal is a clean CLI that wraps stable legacy scripts while staying reliable in short-lived Colab GPU sessions.

Key ideas:
- Colab-first: designed for frequent disconnects and limited storage.
- Terminal-only: no notebooks required.
- Resume-safe MD: OpenMM runs use checkpoint-based chunks.

## Installation (Colab Terminal)

All setup should be done in the Colab Terminal (not notebook cells). These steps also work on any Linux machine.

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
```

Quick sanity check (optional):

```bash
python - << 'PY'
from openmm import Platform
print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PY
```

### 3. Modeller CPU Environment + License

Create a Modeller environment and set your license key:

```bash
conda create -y -n modeller_env -c conda-forge python=3.10 modeller biopython
conda activate modeller_env
```

Get a free academic license at the Modeller site, then update the config:

```bash
LICENSE_KEY="YOUR_LICENSE_KEY"
CONFIG="$HOME/miniforge/envs/modeller_env/lib/modeller-10.8/modlib/modeller/config.py"
sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
python - << 'PY'
import modeller
print('Modeller OK, version', modeller.__version__)
PY
```

## Install ColabMDA (Editable)

From the repository root:

```bash
pip install -e .
```

## Quickstart (Terminal-Only)

All workflows are exposed through the `colabmda` CLI.

### OpenMM (GPU) Workflow

Wild-type (4ldj) example:

```bash
colabmda openmm prep --pdb-id 4ldj
colabmda openmm run --pdb-id 4ldj --total-ns 1 --traj-interval 100 --checkpoint-ps 100
colabmda openmm merge --pdb-id 4ldj
colabmda openmm analysis --pdb-id 4ldj
```

Mutant example (G12C). First generate a mutant PDB (e.g., with Modeller), then run OpenMM using the local file:

```bash
colabmda modeller mutate --pdb-in 4ldj/4ldj_cleaned.pdb --chain A --mut G12C --outdir-mut 4ldj_mut
colabmda openmm prep --pdb-file 4ldj_mut/4ldj_G12C.pdb --name 4ldj_g12c --outdir 4ldj_g12c
colabmda openmm run --workdir 4ldj_g12c --name 4ldj_g12c --total-ns 1 --traj-interval 100 --checkpoint-ps 100
colabmda openmm merge --pdb-dir 4ldj_g12c
colabmda openmm analysis --pdb-dir 4ldj_g12c
```

Use a local PDB file instead of `--pdb-id` when needed:

```bash
colabmda openmm prep --pdb-file /content/drive/MyDrive/4ldj.pdb --name 4ldj --outdir 4ldj
```

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

- This repo preserves legacy research scripts; the CLI wraps them without refactoring.
- Keep dependencies minimal in `pyproject.toml`; OpenMM and Modeller are installed via Conda environments.

If anything in the setup is unclear, open an issue with your exact terminal output.
