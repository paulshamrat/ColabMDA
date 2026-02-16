#!/usr/bin/env bash
set -euo pipefail

# One-command Colab bootstrap using the proven legacy install pattern:
# - Installs Miniforge (if missing)
# - Installs OpenMM stack in conda base environment
# - Installs ColabMDA
# - Verifies GPU/OpenMM/CLI
#
# Usage:
#   bash bootstrap_colab_openmm_gpu.sh
#   bash bootstrap_colab_openmm_gpu.sh <release_tag>
#   WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh <release_tag>
#
# Example:
#   bash bootstrap_colab_openmm_gpu.sh latest

REPO="${COLABMDA_REPO:-paulshamrat/ColabMDA}"
TAG="${1:-latest}"
WITH_MODELLER="${WITH_MODELLER:-0}"
MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"
INSTALL_DIR="${INSTALL_DIR:-/content/colabmda}"
WORK_DIR="${WORK_DIR:-/content/work}"
DRIVE_RUNS_DIR="${DRIVE_RUNS_DIR:-/content/drive/MyDrive/openmm}"
MODELLER_ENV_NAME="${MODELLER_ENV_NAME:-modeller_env}"

echo "[STEP] GPU check"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi || true
else
  echo "[WARN] nvidia-smi not found. Ensure Colab runtime type is GPU."
fi

echo "[STEP] Ensure Miniforge at ${MINIFORGE_DIR}"
if [[ ! -x "${MINIFORGE_DIR}/bin/conda" ]]; then
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p "${MINIFORGE_DIR}"
fi

export PATH="${MINIFORGE_DIR}/bin:${PATH}"
source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
conda activate base

echo "[STEP] Ensure mamba in base"
conda install -y -n base -c conda-forge mamba

echo "[STEP] Install OpenMM stack in conda base"
mamba install -y -n base -c conda-forge cudatoolkit=11.8 openmm openmmtools
conda install -y -n base -c conda-forge pdbfixer || pip install pdbfixer
mamba install -y -n base -c conda-forge mdanalysis mdtraj numpy matplotlib biopython
python -m pip install --upgrade pip

if [[ "${WITH_MODELLER}" != "1" && -t 0 ]]; then
  read -r -p "Install MODELLER environment as well? [y/N]: " _ans
  case "${_ans}" in
    y|Y|yes|YES) WITH_MODELLER="1" ;;
  esac
fi

if [[ "${WITH_MODELLER}" == "1" ]]; then
  echo "[STEP] Install MODELLER environment (${MODELLER_ENV_NAME})"
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
    echo "[WARN] Obtain a key at: https://salilab.org/modeller/"
    echo "[WARN] Set it later: export KEY_MODELLER='YOUR_MODELLER_LICENSE_KEY'"
  else
    conda env config vars set KEY_MODELLER="${KEY_MODELLER}" >/dev/null || true
    if ! grep -q 'KEY_MODELLER=' "$HOME/.bashrc" 2>/dev/null; then
      echo "export KEY_MODELLER='${KEY_MODELLER}'" >> "$HOME/.bashrc"
    fi

    # Patch MODELLER config.py directly as fallback for builds that ignore env-only setup.
    cfg_glob="${MINIFORGE_DIR}/envs/${MODELLER_ENV_NAME}/lib/modeller-*/modlib/modeller/config.py"
    cfg_updated=0
    for cfg in $cfg_glob; do
      if [[ -f "$cfg" ]]; then
        sed -i "s/^license *=.*/license = r'${KEY_MODELLER}'/" "$cfg"
        cfg_updated=1
        echo "[INFO] Updated MODELLER license in: $cfg"
      fi
    done
    if [[ "$cfg_updated" -eq 0 ]]; then
      echo "[WARN] Could not locate MODELLER config.py for direct license patch."
    fi

    echo "[INFO] KEY_MODELLER configured for ${MODELLER_ENV_NAME}"
  fi

  conda activate base
fi

echo "[STEP] Install ColabMDA (${TAG})"
mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

WHEEL_URL="$(
python - "${REPO}" "${TAG}" <<'PY'
import json
import sys
import urllib.request

repo = sys.argv[1]
tag = sys.argv[2]
api = (
    f"https://api.github.com/repos/{repo}/releases/latest"
    if tag == "latest"
    else f"https://api.github.com/repos/{repo}/releases/tags/{tag}"
)

try:
    with urllib.request.urlopen(api, timeout=15) as r:
        data = json.load(r)

    wheels = [
        a["browser_download_url"]
        for a in data.get("assets", [])
        if a.get("name", "").endswith(".whl")
    ]
    print(wheels[0] if wheels else "")
except Exception:
    print("")
PY
)"

if [[ -n "${WHEEL_URL}" ]]; then
  echo "[INFO] Downloading release wheel: ${WHEEL_URL}"
  if curl -fsSL "${WHEEL_URL}" -o colabmda.whl; then
    python -m pip install --upgrade ./colabmda.whl
  else
    echo "[WARN] Release wheel download failed."
  fi
else
  echo "[WARN] Release wheel not found for ${TAG}."
fi

if ! command -v colabmda >/dev/null 2>&1; then
  echo "[INFO] Fallback: trying PyPI package"
  if ! python -m pip install --upgrade colabmda; then
    echo "[WARN] PyPI install failed."
  fi
fi

if ! command -v colabmda >/dev/null 2>&1; then
  if [[ "${TAG}" == "latest" ]]; then
    REF="main"
  else
    REF="${TAG}"
  fi
  echo "[INFO] Fallback: installing from GitHub source (${REF})"
  python -m pip install --upgrade "git+https://github.com/${REPO}.git@${REF}"
fi

echo "[STEP] Validate OpenMM + analysis libs + ColabMDA"
python - <<'PY'
from openmm import Platform
platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print("OpenMM platforms:", platforms)
import MDAnalysis, mdtraj, Bio
print("MDAnalysis:", MDAnalysis.__version__)
print("MDTraj:", mdtraj.__version__)
print("Biopython:", Bio.__version__)
PY

colabmda --help >/dev/null

mkdir -p "${WORK_DIR}" "${DRIVE_RUNS_DIR}"

echo "[OK] Bootstrap complete"
echo "[OK] Env: base"
echo "[OK] Local workdir: ${WORK_DIR}"
echo "[OK] Drive run dir: ${DRIVE_RUNS_DIR}"
if [[ "${WITH_MODELLER}" == "1" ]]; then
  echo "[OK] MODELLER env: ${MODELLER_ENV_NAME}"
fi
echo "[NEXT] Run: conda activate base && cd ${WORK_DIR}"
