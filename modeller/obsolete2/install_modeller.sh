#!/usr/bin/env bash
set -euo pipefail

## CONFIGURATION
ROOT="$HOME/miniforge"
ENV_NAME="modeller_env"
LICENSE_KEY="MODELIRANJE"   # ← replace with your actual key
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
# set license key as an env var so it never prompts
conda env config vars set KEY_MODELLER="$LICENSE_KEY"
mamba install -y modeller

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
