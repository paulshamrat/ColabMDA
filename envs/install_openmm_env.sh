#!/usr/bin/env bash
set -euo pipefail

ENV_NAME=${ENV_NAME:-openmm_gpu}
CUDA_VERSION=${CUDA_VERSION:-11.8}
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install Miniforge/Conda and try again." >&2
  exit 1
fi

# Ensure conda is initialized in this shell
CONDA_BASE=$(conda info --base)
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

if command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
else
  SOLVER="conda"
  conda install -y -n base -c conda-forge conda-libmamba-solver
  conda config --set solver libmamba
fi

$SOLVER create -y -n "$ENV_NAME" -c conda-forge --override-channels \
  python=3.10 \
  openmm \
  openmmtools \
  pdbfixer \
  mdanalysis \
  mdtraj \
  numpy \
  matplotlib \
  biopython \
  "cuda-version=${CUDA_VERSION}"

conda activate "$ENV_NAME"

if [[ -f "$REPO_ROOT/pyproject.toml" ]]; then
  (cd "$REPO_ROOT" && pip install -e .)
else
  echo "WARNING: pyproject.toml not found in $REPO_ROOT. Skipping pip install -e ." >&2
fi

echo "OpenMM environment '$ENV_NAME' is ready."
