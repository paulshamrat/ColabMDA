#!/usr/bin/env bash
set -euo pipefail

## CONFIGURATION
ROOT="$HOME/miniforge"
ENV_NAME="modeller_env"
LICENSE_KEY="${KEY_MODELLER:-${MODELLER_LICENSE:-}}"
LOG="$PWD/install_modeller.log"

exec > >(tee -a "$LOG") 2>&1

echo "===== Conda MODELLER Install Start: $(date) ====="

# 1) Bootstrap Miniforge if needed
if [ ! -x "$ROOT/bin/conda" ]; then
  echo "--- Installing Miniforge"
  wget -qO ~/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash ~/miniforge.sh -b -p "$ROOT"
  rm ~/miniforge.sh
fi
export PATH="$ROOT/bin:$PATH"

# 2) Initialize conda in this shell
eval "$(conda shell.bash hook)"

# 3) Create or update the modeller_env
if conda env list | grep -q "${ENV_NAME}"; then
  echo "--- Updating existing environment: $ENV_NAME"
  conda activate "$ENV_NAME"
else
  echo "--- Creating environment: $ENV_NAME"
  conda create -y -n "$ENV_NAME" python=3.10
  conda activate "$ENV_NAME"
fi

# 4) Configure channels & install MODELLER
echo "--- Configuring channels and installing MODELLER"
conda config --env --add channels salilab
conda config --env --add channels conda-forge

if [[ -z "${LICENSE_KEY}" ]]; then
  if [[ -t 0 ]]; then
    read -r -p "Enter MODELLER license key (KEY_MODELLER): " LICENSE_KEY
  else
    echo "❌ KEY_MODELLER is not set and no interactive prompt is available"
    echo "Set it before running, e.g.: export KEY_MODELLER='MODELIRANJE'"
    exit 1
  fi
fi

if [[ -z "${LICENSE_KEY}" ]]; then
  echo "❌ Empty MODELLER license key."
  exit 1
fi

# Set license key as an env var so MODELLER does not prompt later.
conda env config vars set KEY_MODELLER="$LICENSE_KEY"
if ! grep -q 'KEY_MODELLER=' "$HOME/.bashrc" 2>/dev/null; then
  echo "export KEY_MODELLER='${LICENSE_KEY}'" >> "$HOME/.bashrc"
fi

mamba install -y modeller

# Patch MODELLER config.py directly; required on some builds.
cfg_glob="${ROOT}/envs/${ENV_NAME}/lib/modeller-*/modlib/modeller/config.py"
cfg_updated=0
for cfg in $cfg_glob; do
  if [[ -f "$cfg" ]]; then
    sed -i "s/^license *=.*/license = r'${LICENSE_KEY}'/" "$cfg"
    cfg_updated=1
    echo "--- Updated MODELLER license in: $cfg"
  fi
done
if [[ "$cfg_updated" -eq 0 ]]; then
  echo "WARNING: Could not locate MODELLER config.py for direct license patch."
fi

# 5) Verify
echo "--- Verifying installation"
if command -v mod10.7 >/dev/null 2>&1; then
  echo "✅ mod10.7 is available at $(which mod10.7)"
  mod10.7 | head -1
else
  echo "❌ mod10.7 not found on PATH"
  exit 1
fi

echo "===== Installation Complete: $(date) ====="
