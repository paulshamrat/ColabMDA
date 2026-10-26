#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# ColabMDA Universal Bootstrap Script (Google Colab & Local Workstation)
# ==============================================================================
#
# Usage:
#   bash scripts/setup.sh
#   bash scripts/setup.sh <release_tag>
#   WITH_MODELLER=1 bash scripts/setup.sh
#   WITH_MODELLER=1 KEY_MODELLER="YOUR_LICENSE_KEY" bash scripts/setup.sh
#
# What this script does:
#   1. Detects environment automatically (Google Colab VM vs Local Workstation)
#   2. Verifies GPU availability (Tesla T4 / A100 / L4 / local NVIDIA or CPU fallback)
#   3. Installs Miniforge conda manager in $HOME/miniforge3 (if missing)
#   4. Provisions core simulation environment (OpenMM, AmberTools, MDAnalysis)
#   5. Configures MODELLER (enabled by default locally, unattended on Colab)
#   6. Installs PyMOL & google-colab-cli on local workstation environments
#   7. Installs ColabMDA and registers CLI entrypoints (colabmda and cmda)
#   8. Runs an automated post-install health check
# ==============================================================================

REPO="${COLABMDA_REPO:-paulshamrat/ColabMDA}"
TAG="${1:-latest}"
OPENMM_CUDA_VERSION="${COLABMDA_CUDA_VERSION:-13.0}"
ENV_NAME="${ENV_NAME:-cmda}"
MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"

# Detect execution platform
IS_COLAB=0
if [[ -d "/content" ]]; then
  IS_COLAB=1
fi

if [[ "${IS_COLAB}" == "1" ]]; then
  TARGET_ENV="Google Colab (Cloud GPU)"
  WITH_MODELLER="${WITH_MODELLER:-0}"
  WITH_PYMOL="${WITH_PYMOL:-0}"
  WITH_COLAB_CLI="${WITH_COLAB_CLI:-0}"
else
  TARGET_ENV="Local Workstation / Linux"
  WITH_MODELLER="${WITH_MODELLER:-1}"
  WITH_PYMOL="${WITH_PYMOL:-1}"
  WITH_COLAB_CLI="${WITH_COLAB_CLI:-1}"
fi

WITH_LIGAND="${WITH_LIGAND:-0}"
WITH_AMBERTOOLS="${WITH_AMBERTOOLS:-1}"

DRIVE_RUNS_DIR="${DRIVE_RUNS_DIR:-}"
if [[ -z "${DRIVE_RUNS_DIR}" ]] && [[ -d "/content/drive/MyDrive" ]]; then
  DRIVE_RUNS_DIR="/content/drive/MyDrive/research/projects"
fi

SCRIPT_DIR=""
if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd 2>/dev/null || true)"
fi
REPO_DIR=""
if [[ -n "${SCRIPT_DIR}" ]]; then
  REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd 2>/dev/null || true)"
fi
if [[ -z "${COLABMDA_SOURCE:-}" ]] && [[ -n "${REPO_DIR}" ]] && [[ -f "${REPO_DIR}/pyproject.toml" ]]; then
  COLABMDA_SOURCE="${REPO_DIR}"
fi

INSTALL_REF="${INSTALL_REF:-}"
if [[ -z "${INSTALL_REF}" ]]; then
  if [[ "${TAG}" == "latest" ]]; then
    INSTALL_REF="dev1"
  else
    INSTALL_REF="${TAG}"
  fi
fi

echo "=============================================================================="
echo " Starting ColabMDA Environment Bootstrap (${TAG})"
echo " Target Platform:    ${TARGET_ENV}"
echo " Environment Name:   ${ENV_NAME}"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 1. GPU Check
# ------------------------------------------------------------------------------
echo "[STEP 1/6] Checking GPU availability..."
HAS_GPU=0
if command -v nvidia-smi >/dev/null 2>&1; then
  HAS_GPU=1
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader || true
else
  if [[ "${IS_COLAB}" == "1" ]]; then
    echo "[WARN] nvidia-smi not detected on Colab. Ensure runtime type is set to GPU (T4/A100)."
  else
    echo "[INFO] No NVIDIA GPU detected locally. Configuring for CPU/OpenCL OpenMM..."
  fi
fi

# ------------------------------------------------------------------------------
# 2. Miniforge & Mamba Installation
# ------------------------------------------------------------------------------
echo "[STEP 2/6] Ensuring Miniforge at ${MINIFORGE_DIR}..."
if ! command -v conda >/dev/null 2>&1 && [[ ! -x "${MINIFORGE_DIR}/bin/conda" ]]; then
  echo "[INFO] Downloading and installing Miniforge..."
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p "${MINIFORGE_DIR}"
  rm -f /tmp/miniforge.sh
fi

if [[ -f "${MINIFORGE_DIR}/etc/profile.d/conda.sh" ]]; then
  export PATH="${MINIFORGE_DIR}/bin:${PATH}"
  # shellcheck source=/dev/null
  source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
fi

conda activate base

if ! command -v mamba >/dev/null 2>&1; then
  echo "[INFO] Installing mamba in base environment..."
  conda install -y -n base -c conda-forge mamba
fi

# ------------------------------------------------------------------------------
# 3. Core OpenMM + AmberTools Simulation Environment
# ------------------------------------------------------------------------------
echo "[STEP 3/6] Setting up simulation environment: ${ENV_NAME}..."
SPEC=$(mktemp /tmp/colabmda-environment-XXXXXX.yml)
trap 'rm -f "$SPEC"' EXIT
if [[ -n "${COLABMDA_SOURCE:-}" && -f "${COLABMDA_SOURCE}/environment.yml" ]]; then
  cp "$COLABMDA_SOURCE/environment.yml" "$SPEC"
elif [[ -n "${SCRIPT_DIR:-}" && -f "${SCRIPT_DIR}/../environment.yml" ]]; then
  cp "${SCRIPT_DIR}/../environment.yml" "$SPEC"
else
  wget -q "https://raw.githubusercontent.com/${REPO}/${INSTALL_REF}/environment.yml" -O "$SPEC"
fi

if [[ "${HAS_GPU}" == "1" ]]; then
  [[ "$OPENMM_CUDA_VERSION" =~ ^[0-9]+\.[0-9]+$ ]] || { echo 'Invalid CUDA version' >&2; exit 2; }
  sed -i -E "s/cuda-version=[0-9.]+/cuda-version=${OPENMM_CUDA_VERSION}/" "$SPEC"
else
  # Strip cuda-version requirement for non-GPU machines
  sed -i '/cuda-version=/d' "$SPEC"
fi

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
# 4. MODELLER, PyMOL, and Google Colab CLI Setup
# ------------------------------------------------------------------------------
if [[ "${WITH_MODELLER}" != "1" && -n "${KEY_MODELLER:-}" ]]; then
  WITH_MODELLER="1"
fi

if [[ "${WITH_MODELLER}" == "1" ]]; then
  echo "[STEP 4/6] Setting up MODELLER..."
  mamba install -y -q -n "${ENV_NAME}" -c salilab -c conda-forge modeller
  conda activate "${ENV_NAME}"

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
else
  echo "[STEP 4/6] Skipping MODELLER (pass WITH_MODELLER=1 to enable)."
fi

if [[ "${WITH_PYMOL}" == "1" ]]; then
  echo "[INFO] Installing Open-Source PyMOL for molecular visualization..."
  mamba install -y -q -n "${ENV_NAME}" -c conda-forge pymol-open-source || true
fi

if [[ "${WITH_COLAB_CLI}" == "1" ]]; then
  echo "[INFO] Installing google-colab-cli for cloud coordination..."
  python -m pip install --quiet google-colab-cli || true
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

# Register global CLI shortcuts in /usr/local/bin if writable
if [[ -w "/usr/local/bin" ]]; then
  ln -sf "${MINIFORGE_DIR}/envs/${ENV_NAME}/bin/cmda" /usr/local/bin/cmda 2>/dev/null || true
  ln -sf "${MINIFORGE_DIR}/envs/${ENV_NAME}/bin/colabmda" /usr/local/bin/colabmda 2>/dev/null || true
fi

# Auto-activate in ~/.bashrc so new subshells have cmda immediately available
if ! grep -q "conda activate ${ENV_NAME}" "${HOME}/.bashrc" 2>/dev/null; then
  {
    echo ""
    echo "# ColabMDA environment auto-activation"
    echo "if [[ -f \"${MINIFORGE_DIR}/etc/profile.d/conda.sh\" ]]; then"
    echo "  source \"${MINIFORGE_DIR}/etc/profile.d/conda.sh\""
    echo "  conda activate ${ENV_NAME}"
    echo "fi"
  } >> "${HOME}/.bashrc"
fi

# ------------------------------------------------------------------------------
# 6. Automated Post-Install Health Check (Self-Test)
# ------------------------------------------------------------------------------
echo "[STEP 6/6] Running automated post-install health check..."
python - <<'PY'
import sys
from openmm import Platform
from openmm.app import ForceField
import colabmda

# 1. Platform validation
platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print(f"[HEALTH] Available OpenMM platforms: {platforms}")
if "CUDA" in platforms:
    print("[HEALTH] CUDA platform: OK")
else:
    print("[HEALTH] OpenMM CPU/Reference platform: OK")

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

if command -v tleap >/dev/null 2>&1 || [[ -x "${MINIFORGE_DIR}/envs/${ENV_NAME}/bin/tleap" ]]; then
  echo "[HEALTH] AmberTools tleap: OK"
fi

if python -c "import modeller" >/dev/null 2>&1; then
  echo "[HEALTH] MODELLER: OK"
fi

if python -c "import pymol" >/dev/null 2>&1; then
  echo "[HEALTH] PyMOL: OK"
fi

if command -v colab >/dev/null 2>&1 || [[ -x "${MINIFORGE_DIR}/envs/${ENV_NAME}/bin/colab" ]]; then
  echo "[HEALTH] Google Colab CLI: OK"
fi

if [[ -n "$DRIVE_RUNS_DIR" ]]; then
  mkdir -p "$DRIVE_RUNS_DIR"
fi

echo "=============================================================================="
echo " [OK] ColabMDA Environment Successfully Installed!"
echo " Target Platform:    ${TARGET_ENV}"
echo " Env Name:           ${ENV_NAME}"
if [[ -n "${DRIVE_RUNS_DIR}" ]]; then
  echo " Google Drive:       ${DRIVE_RUNS_DIR}"
fi
echo ""
echo " IMPORTANT: Activate the environment in this terminal session:"
echo "   source ${MINIFORGE_DIR}/etc/profile.d/conda.sh"
echo "   conda activate ${ENV_NAME}"
echo "=============================================================================="
