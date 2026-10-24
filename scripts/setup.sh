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
INSTALL_DIR="${INSTALL_DIR:-/content/colabmda}"
WORK_DIR="${WORK_DIR:-/content/work}"
DRIVE_RUNS_DIR="${DRIVE_RUNS_DIR:-}"
OPENMM_CUDA_VERSION="${COLABMDA_CUDA_VERSION:-12.4}"
ENV_NAME="${ENV_NAME:-openmm_env}"
MODELLER_ENV_NAME="${MODELLER_ENV_NAME:-modeller_env}"

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
PACKAGES=(
  python=3.10
  "cuda-version=${OPENMM_CUDA_VERSION}"
  "openmm>=8.2"
  pdbfixer
  pdb2pqr
  propka
  mdanalysis
  mdtraj
  numpy
  matplotlib
  biopython
)

if [[ "${WITH_AMBERTOOLS}" == "1" ]]; then
  PACKAGES+=(ambertools)
fi

if [[ "${WITH_LIGAND}" == "1" ]]; then
  PACKAGES+=(openff-toolkit openmmforcefields)
fi

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[INFO] Updating existing environment: ${ENV_NAME}..."
  mamba install -y -n "${ENV_NAME}" -c conda-forge "${PACKAGES[@]}"
else
  echo "[INFO] Creating new environment: ${ENV_NAME}..."
  mamba create -y -n "${ENV_NAME}" -c conda-forge "${PACKAGES[@]}"
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
  mamba create -y -n "${MODELLER_ENV_NAME}" -c salilab -c conda-forge python=3.10 modeller biopython
  conda activate "${MODELLER_ENV_NAME}"
  python -m pip install --upgrade pip

  if [[ -z "${KEY_MODELLER:-}" && -n "${MODELLER_LICENSE:-}" ]]; then
    export KEY_MODELLER="${MODELLER_LICENSE}"
  fi
  if [[ -z "${KEY_MODELLER:-}" && -t 0 ]]; then
    read -r -p "Enter MODELLER license key (KEY_MODELLER): " KEY_MODELLER
    export KEY_MODELLER
  fi
  if [[ -z "${KEY_MODELLER:-}" ]]; then
    echo "[WARN] MODELLER installed but KEY_MODELLER not set."
    echo "[WARN] Obtain a free academic key at: https://salilab.org/modeller/"
    echo "[WARN] Set it later: export KEY_MODELLER='YOUR_LICENSE_KEY'"
  else
    conda env config vars set KEY_MODELLER="${KEY_MODELLER}" >/dev/null || true
    if ! grep -q 'KEY_MODELLER=' "$HOME/.bashrc" 2>/dev/null; then
      echo "export KEY_MODELLER='${KEY_MODELLER}'" >> "$HOME/.bashrc"
    fi
    # Patch MODELLER config.py directly as fallback
    for cfg in ${MINIFORGE_DIR}/envs/${MODELLER_ENV_NAME}/lib/modeller-*/modlib/modeller/config.py; do
      if [[ -f "$cfg" ]]; then
        sed -i "s/^license *=.*/license = r'${KEY_MODELLER}'/" "$cfg"
        echo "[INFO] Injected license into: $cfg"
      fi
    done
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
mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

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
print(f"[HEALTH] ColabMDA v{colabmda.__version__}: OK")
PY

if command -v tleap >/dev/null 2>&1; then
  echo "[HEALTH] AmberTools tleap: $(which tleap)"
fi

mkdir -p "${WORK_DIR}"
[[ -z "$DRIVE_RUNS_DIR" ]] || mkdir -p "$DRIVE_RUNS_DIR"
true

echo "=============================================================================="
echo " [OK] ColabMDA GPU Environment Successfully Installed!"
echo " Env Name:       ${ENV_NAME}"
echo " Working Dir:    ${WORK_DIR}"
echo " Google Drive:   ${DRIVE_RUNS_DIR}"
echo " Activate via:   conda activate ${ENV_NAME} && cd ${WORK_DIR}"
echo "=============================================================================="
