#!/usr/bin/env bash
set -euo pipefail

OPENMMDL_REPO_URL="${OPENMMDL_REPO_URL:-https://github.com/wolberlab/OpenMMDL.git}"
OPENMMDL_DIR="/content/OpenMMDL"
MAMBA_ROOT_PREFIX="/content/micromamba"

echo "[1/7] Clone/Open OpenMMDL at ${OPENMMDL_DIR}"
if [[ -d "${OPENMMDL_DIR}/.git" ]]; then
  echo "OpenMMDL already exists, pulling latest changes."
  git -C "${OPENMMDL_DIR}" pull --ff-only
else
  git clone "${OPENMMDL_REPO_URL}" "${OPENMMDL_DIR}"
fi

echo "[2/7] Install micromamba"
cd /content
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /tmp bin/micromamba
mkdir -p "${MAMBA_ROOT_PREFIX}/bin"
cp /tmp/bin/micromamba "${MAMBA_ROOT_PREFIX}/bin/"

echo "[3/7] Activate micromamba shell"
export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"

echo "[4/7] Recreate env: openmmdl"
if micromamba env list | awk '{print $1}' | grep -qx "openmmdl"; then
  micromamba env remove -y -n openmmdl
fi
micromamba env create -y -n openmmdl -f "${OPENMMDL_DIR}/environment.yml"

echo "[5/7] Install OpenMMDL package"
micromamba activate openmmdl
cd "${OPENMMDL_DIR}"
pip install .

echo "[6/7] Basic checks"
nvidia-smi || true
openmmdl --help >/tmp/openmmdl_help.txt
head -n 20 /tmp/openmmdl_help.txt

echo "[7/7] OpenMM CUDA check"
python - <<'PY'
import openmm as mm
platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print("OpenMM version:", mm.__version__)
print("Platforms:", platforms)
print("CUDA available:", "CUDA" in platforms)
PY

echo "Install completed."
