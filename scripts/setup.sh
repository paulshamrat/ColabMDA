#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# ColabMDA Unified Google Colab GPU Bootstrap Script
# ==============================================================================
#
# Usage:
#   bash scripts/setup.sh
#   bash scripts/setup.sh <release_tag>
#   WITH_MODELLER=1 bash scripts/setup.sh
#   WITH_MODELLER=1 KEY_MODELLER="YOUR_LICENSE_KEY" bash scripts/setup.sh
#
# What this script does:
#   1. Checks NVIDIA GPU presence (Tesla T4 / A100 / L4)
#   2. Installs Miniforge conda manager in $HOME/miniforge3 (if missing)
#   3. Provisions the core MD environment (OpenMM + CUDA, AmberTools, MDAnalysis)
#   4. Configures MODELLER (if WITH_MODELLER=1 or KEY_MODELLER is provided)
#   5. Installs ColabMDA and registers CLI entrypoints (colabmda and cmda)
#   6. Runs an automated GPU self-test verifying CUDA platform & force fields
# ==============================================================================

REPO="${COLABMDA_REPO:-paulshamrat/ColabMDA}"
TAG="${1:-latest}"
WITH_MODELLER="${WITH_MODELLER:-0}"
WITH_LIGAND="${WITH_LIGAND:-0}"
WITH_AMBERTOOLS="${WITH_AMBERTOOLS:-1}"
MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"
DRIVE_RUNS_DIR="${DRIVE_RUNS_DIR:-}"
if [[ -z "${DRIVE_RUNS_DIR}" ]] && [[ -d "/content/drive/MyDrive" ]]; then
  DRIVE_RUNS_DIR="/content/drive/MyDrive/research/projects"
fi
OPENMM_CUDA_VERSION="${COLABMDA_CUDA_VERSION:-13.0}"
ENV_NAME="${ENV_NAME:-colabmda}"
MODELLER_ENV_NAME="${MODELLER_ENV_NAME:-modeller_env}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [[ -z "${COLABMDA_SOURCE:-}" ]] && [[ -f "${REPO_DIR}/pyproject.toml" ]]; then
  COLABMDA_SOURCE="${REPO_DIR}"
fi

INSTALL_REF="${INSTALL_REF:-}"
if [[ -z "${INSTALL_REF}" ]]; then
  if [[ "${TAG}" == "latest" ]]; then
    INSTALL_REF="main"
  else
    INSTALL_REF="${TAG}"
  fi
fi

echo "=============================================================================="
echo " Starting ColabMDA GPU Environment Bootstrap (${TAG})"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 1. GPU Check
# ------------------------------------------------------------------------------
echo "[STEP 1/6] Checking GPU availability..."
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader || true
else
  echo "[WARN] nvidia-smi not detected. Ensure your Colab runtime is set to GPU (T4/A100)."
fi

# ------------------------------------------------------------------------------
# 2. Miniforge & Mamba Installation
# ------------------------------------------------------------------------------
echo "[STEP 2/6] Ensuring Miniforge at ${MINIFORGE_DIR}..."
if [[ ! -x "${MINIFORGE_DIR}/bin/conda" ]]; then
  echo "[INFO] Downloading and installing Miniforge..."
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p "${MINIFORGE_DIR}"
  rm -f /tmp/miniforge.sh
fi

export PATH="${MINIFORGE_DIR}/bin:${PATH}"
# shellcheck source=/dev/null
source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
conda activate base

if ! command -v mamba >/dev/null 2>&1; then
  echo "[INFO] Installing mamba in base environment..."
  conda install -y -n base -c conda-forge mamba
fi

# ------------------------------------------------------------------------------
# 3. Core OpenMM + AmberTools Simulation Environment
# ------------------------------------------------------------------------------
echo "[STEP 3/6] Setting up simulation environment: ${ENV_NAME}..."
# The readable environment specification is the single core dependency list.
SPEC=$(mktemp /tmp/colabmda-environment-XXXXXX.yml)
trap 'rm -f "$SPEC"' EXIT
if [[ -f "${COLABMDA_SOURCE:-}/environment.yml" ]]; then
  cp "$COLABMDA_SOURCE/environment.yml" "$SPEC"
elif [[ -f "$(dirname -- "${BASH_SOURCE[0]}")/../environment.yml" ]]; then
  cp "$(dirname -- "${BASH_SOURCE[0]}")/../environment.yml" "$SPEC"
else
  wget -q "https://raw.githubusercontent.com/${REPO}/${INSTALL_REF}/environment.yml" -O "$SPEC"
fi
[[ "$OPENMM_CUDA_VERSION" =~ ^[0-9]+\.[0-9]+$ ]] || { echo 'Invalid CUDA version' >&2; exit 2; }
sed -i -E "s/cuda-version=[0-9.]+/cuda-version=${OPENMM_CUDA_VERSION}/" "$SPEC"
if [[ "$WITH_AMBERTOOLS" != 1 ]]; then sed -i '/^  - ambertools$/d' "$SPEC"; fi
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  mamba env update -y -q -n "$ENV_NAME" -f "$SPEC"
else
  mamba env create -y -q -n "$ENV_NAME" -f "$SPEC"
fi
if [[ "$WITH_LIGAND" == 1 ]]; then
  mamba install -y -q -n "$ENV_NAME" -c conda-forge openff-toolkit openmmforcefields
fi

conda activate "${ENV_NAME}"
python -m pip install --upgrade pip

# ------------------------------------------------------------------------------
# 4. Optional MODELLER Setup
# ------------------------------------------------------------------------------
if [[ "${WITH_MODELLER}" != "1" && -n "${KEY_MODELLER:-}" ]]; then
  WITH_MODELLER="1"
fi

if [[ "${WITH_MODELLER}" != "1" && -t 0 ]]; then
  read -r -p "Install MODELLER environment as well? [y/N]: " _ans
  case "${_ans}" in
    y|Y|yes|YES) WITH_MODELLER="1" ;;
  esac
fi

if [[ "${WITH_MODELLER}" == "1" ]]; then
  echo "[STEP 4/6] Setting up MODELLER environment (${MODELLER_ENV_NAME})..."
  conda config --add channels salilab || true
  if conda env list | awk '{print $1}' | grep -qx "$MODELLER_ENV_NAME"; then
    mamba install -y -q -n "$MODELLER_ENV_NAME" -c salilab -c conda-forge python=3.10 modeller biopython
  else
    mamba create -y -q -n "$MODELLER_ENV_NAME" -c salilab -c conda-forge python=3.10 modeller biopython
  fi
  conda activate "${MODELLER_ENV_NAME}"
  python -m pip install --upgrade pip

  if [[ -z "${KEY_MODELLER:-}" && -n "${MODELLER_LICENSE:-}" ]]; then
    export KEY_MODELLER="${MODELLER_LICENSE}"
  fi
  if [[ -z "${KEY_MODELLER:-}" && -t 0 ]]; then
    read -r -s -p "Enter MODELLER license key (KEY_MODELLER): " KEY_MODELLER
    printf '\n'
    export KEY_MODELLER
  fi
  if [[ -z "${KEY_MODELLER:-}" ]]; then
    echo "[WARN] MODELLER installed but KEY_MODELLER not set."
    echo "[WARN] Obtain a free academic key at: https://salilab.org/modeller/"
    echo "[WARN] Set it later: export KEY_MODELLER='YOUR_LICENSE_KEY'"
  else
    # Read the key from the process environment, never a command-line argument.
    python - <<'PY_LICENSE'
import json
import os
import re
import tempfile
from pathlib import Path

prefix = Path(os.environ["CONDA_PREFIX"])
prefix.chmod(0o700)
key = os.environ["KEY_MODELLER"]

def private_write(path, text):
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(text)
    temporary.chmod(0o600)
    temporary.replace(path)

state = prefix / "conda-meta/state"
settings = json.loads(state.read_text()) if state.exists() else {}
settings.setdefault("env_vars", {})["KEY_MODELLER"] = key
private_write(state, json.dumps(settings, indent=2) + "\n")
for cfg in prefix.glob("lib/modeller-*/modlib/modeller/config.py"):
    text = re.sub(r"^license\s*=.*$", lambda _: "license = " + repr(key), cfg.read_text(), flags=re.M)
    private_write(cfg, text)
print("[INFO] Private MODELLER environment and license files configured.")
PY_LICENSE
    echo "[INFO] KEY_MODELLER successfully configured."
  fi

  conda activate "${ENV_NAME}"
else
  echo "[STEP 4/6] Skipping MODELLER (pass WITH_MODELLER=1 to enable)."
fi

# ------------------------------------------------------------------------------
# 5. ColabMDA Package Installation
# ------------------------------------------------------------------------------
echo "[STEP 5/6] Installing ColabMDA (${TAG})..."

if [[ -n "${COLABMDA_SOURCE:-}" ]]; then
  python -m pip install -e "$COLABMDA_SOURCE"
else
  python -m pip install --upgrade "git+https://github.com/${REPO}.git@${INSTALL_REF}"
fi

# ------------------------------------------------------------------------------
# 6. Automated Post-Install Health Check (Self-Test)
# ------------------------------------------------------------------------------
echo "[STEP 6/6] Running automated post-install GPU health check..."
python - <<'PY'
import sys
from openmm import Platform
from openmm.app import ForceField
import colabmda

# 1. Platform validation
platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print(f"[HEALTH] Available OpenMM platforms: {platforms}")
if "CUDA" not in platforms:
    print("[WARN] CUDA platform not detected in OpenMM! Verify NVIDIA driver.")
else:
    print("[HEALTH] CUDA platform: OK")

# 2. Forcefield validation
try:
    ForceField("amber19-all.xml", "amber19/opc.xml")
    print("[HEALTH] Amber19/OPC ForceField XML: OK")
except Exception as e:
    print(f"[WARN] ForceField XML check warning: {e}")

# 3. Core libraries
import MDAnalysis
print(f"[HEALTH] MDAnalysis v{MDAnalysis.__version__}: OK")
version = getattr(colabmda, "__version__", getattr(colabmda, "VERSION", "0.2.0"))
print(f"[HEALTH] ColabMDA v{version}: OK")
PY

if command -v tleap >/dev/null 2>&1; then
  echo "[HEALTH] AmberTools tleap: $(which tleap)"
fi

if [[ -n "$DRIVE_RUNS_DIR" ]]; then
  mkdir -p "$DRIVE_RUNS_DIR"
fi

echo "=============================================================================="
echo " [OK] ColabMDA GPU Environment Successfully Installed!"
echo " Env Name:       ${ENV_NAME}"
if [[ -n "${DRIVE_RUNS_DIR}" ]]; then
  echo " Google Drive:   ${DRIVE_RUNS_DIR}"
fi
echo " Activate via:   conda activate ${ENV_NAME}"
echo "=============================================================================="
