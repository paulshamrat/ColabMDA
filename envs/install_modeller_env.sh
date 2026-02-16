#!/usr/bin/env bash
set -euo pipefail

ENV_NAME=${ENV_NAME:-modeller_env}
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  if [[ -x "$HOME/miniforge3/bin/conda" ]]; then
    export PATH="$HOME/miniforge3/bin:$PATH"
  fi
fi

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

$SOLVER create -y -n "$ENV_NAME" -c salilab -c conda-forge \
  python=3.10 \
  modeller \
  biopython

conda activate "$ENV_NAME"

if [[ -f "$REPO_ROOT/pyproject.toml" ]]; then
  (cd "$REPO_ROOT" && pip install -e .)
else
  echo "WARNING: pyproject.toml not found in $REPO_ROOT. Skipping pip install -e ." >&2
fi

# License handling:
# Priority: KEY_MODELLER > MODELLER_LICENSE > interactive prompt
if [[ -z "${KEY_MODELLER:-}" && -n "${MODELLER_LICENSE:-}" ]]; then
  export KEY_MODELLER="${MODELLER_LICENSE}"
fi

if [[ -z "${KEY_MODELLER:-}" ]]; then
  if [[ -t 0 ]]; then
    read -r -p "Enter MODELLER license key (KEY_MODELLER): " KEY_MODELLER
    export KEY_MODELLER
  else
    echo "ERROR: KEY_MODELLER is not set and no interactive prompt is available." >&2
    echo "Set it before running, e.g.:" >&2
    echo "  export KEY_MODELLER='MODELIRANJE'" >&2
    exit 1
  fi
fi

if [[ -z "${KEY_MODELLER:-}" ]]; then
  echo "ERROR: empty MODELLER license key." >&2
  exit 1
fi

# Persist key for this conda env and future shells.
conda env config vars set KEY_MODELLER="${KEY_MODELLER}" >/dev/null
if ! grep -q 'KEY_MODELLER=' "$HOME/.bashrc" 2>/dev/null; then
  echo "export KEY_MODELLER='${KEY_MODELLER}'" >> "$HOME/.bashrc"
fi

# Also patch MODELLER config.py directly; some builds ignore env var-only setup.
cfg_glob="${CONDA_PREFIX}/lib/modeller-*/modlib/modeller/config.py"
cfg_updated=0
for cfg in $cfg_glob; do
  if [[ -f "$cfg" ]]; then
    sed -i "s/^license *=.*/license = r'${KEY_MODELLER}'/" "$cfg"
    cfg_updated=1
    echo "Updated MODELLER license in: $cfg"
  fi
done
if [[ "$cfg_updated" -eq 0 ]]; then
  echo "WARNING: Could not locate MODELLER config.py for direct license patch." >&2
fi

echo "KEY_MODELLER configured for env '$ENV_NAME'."
echo "Modeller environment '$ENV_NAME' is ready."
