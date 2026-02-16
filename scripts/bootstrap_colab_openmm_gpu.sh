#!/usr/bin/env bash
set -euo pipefail

# One-command Colab bootstrap:
# - Installs Miniforge (if missing)
# - Creates/updates OpenMM GPU environment
# - Installs ColabMDA from GitHub Release wheel (no repo clone)
# - Verifies GPU/OpenMM/CLI
#
# Usage:
#   bash bootstrap_colab_openmm_gpu.sh
#   bash bootstrap_colab_openmm_gpu.sh <release_tag> <env_name>
#
# Example:
#   bash bootstrap_colab_openmm_gpu.sh latest openmm_gpu

REPO="${COLABMDA_REPO:-paulshamrat/ColabMDA}"
TAG="${1:-latest}"
ENV_NAME="${2:-openmm_gpu}"
MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"
INSTALL_DIR="${INSTALL_DIR:-/content/colabmda}"
WORK_DIR="${WORK_DIR:-/content/work}"
DRIVE_RUNS_DIR="${DRIVE_RUNS_DIR:-/content/drive/MyDrive/openmm_runs}"

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

echo "[STEP] Ensure mamba in base"
conda install -y -n base -c conda-forge mamba

echo "[STEP] Create/update env: ${ENV_NAME}"
mamba create -y -n "${ENV_NAME}" -c conda-forge \
  python=3.10 \
  cudatoolkit=11.8 \
  openmm \
  openmmtools \
  pdbfixer \
  mdanalysis \
  mdtraj \
  numpy \
  matplotlib \
  biopython

conda activate "${ENV_NAME}"
python -m pip install --upgrade pip

echo "[STEP] Install ColabMDA from GitHub Release (${TAG})"
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

with urllib.request.urlopen(api) as r:
    data = json.load(r)

wheels = [
    a["browser_download_url"]
    for a in data.get("assets", [])
    if a.get("name", "").endswith(".whl")
]
if not wheels:
    raise SystemExit("ERROR: no wheel asset found in release.")
print(wheels[0])
PY
)"

curl -fL "${WHEEL_URL}" -o colabmda.whl
python -m pip install --upgrade ./colabmda.whl

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
echo "[OK] Env: ${ENV_NAME}"
echo "[OK] Local workdir: ${WORK_DIR}"
echo "[OK] Drive run dir: ${DRIVE_RUNS_DIR}"
echo "[NEXT] Run: conda activate ${ENV_NAME} && cd ${WORK_DIR}"
