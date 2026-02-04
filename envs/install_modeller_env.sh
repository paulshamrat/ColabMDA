#!/usr/bin/env bash
set -euo pipefail

ENV_NAME=${ENV_NAME:-modeller_env}
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

conda config --add channels salilab

$SOLVER create -y -n "$ENV_NAME" -c conda-forge --override-channels \
  python=3.10 \
  modeller \
  biopython

conda activate "$ENV_NAME"

if [[ -f "$REPO_ROOT/pyproject.toml" ]]; then
  (cd "$REPO_ROOT" && pip install -e .)
else
  echo "WARNING: pyproject.toml not found in $REPO_ROOT. Skipping pip install -e ." >&2
fi

if [[ -z "${KEY_MODELLER:-}" ]]; then
  export KEY_MODELLER="MODELIRANJE"
  echo "NOTE: KEY_MODELLER was not set; defaulting to MODELIRANJE." >&2
else
  echo "KEY_MODELLER is set." >&2
fi

echo "Modeller environment '$ENV_NAME' is ready."
