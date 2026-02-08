# OpenMMDL Colab T4 (from ColabMDA)

This folder is designed to be pushed from your **ColabMDA** repository.
It gives a script-first Google Colab Terminal flow to:

1. Install OpenMMDL with micromamba on Colab T4
2. Run a short protein-ligand simulation smoke test

## Files

- `01_install_openmmdl_colab_t4.sh`
- `02_run_quick_protein_ligand_simulation.sh`
- `openmmdl_protein_ligand_quickrun.py`

## Colab Terminal Steps

Set Colab runtime to **GPU** first.

```bash
cd /content
git clone https://github.com/paulshamrat/ColabMDA.git ColabMDA
cd ColabMDA/openmm/openmmdl_colab_t4
bash 01_install_openmmdl_colab_t4.sh
bash 02_run_quick_protein_ligand_simulation.sh
```

## What gets installed

- micromamba at `/content/micromamba`
- OpenMMDL source at `/content/OpenMMDL`
- conda env: `openmmdl`

## Quick simulation details

- Uses OpenMMDL tutorial input:
  - `5wyz-moe-processed_openMMDL.pdb`
  - `5VF.sdf`
- Runs with `openmmdl simulation ...`
- Uses a short step count for smoke testing on T4
- Output folder: `/content/openmmdl_t4_quickrun`

## Expected checks

- `nvidia-smi` works
- OpenMM sees CUDA platform
- `openmmdl --help` works
- Simulation writes trajectory/log/output files in `/content/openmmdl_t4_quickrun`
