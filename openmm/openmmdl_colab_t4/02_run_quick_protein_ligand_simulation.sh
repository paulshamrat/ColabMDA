#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OPENMMDL_DIR="/content/OpenMMDL"
MAMBA_ROOT_PREFIX="/content/micromamba"
RUN_DIR="/content/openmmdl_t4_quickrun"
TUTORIAL_DIR="${OPENMMDL_DIR}/openmmdl/openmmdl_simulation/tutorial_systems/pdb_path/5wyz_solvent"

if [[ ! -d "${OPENMMDL_DIR}" ]]; then
  echo "ERROR: ${OPENMMDL_DIR} not found. Run 01_install_openmmdl_colab_t4.sh first."
  exit 1
fi

if [[ ! -f "${SCRIPT_DIR}/openmmdl_protein_ligand_quickrun.py" ]]; then
  echo "ERROR: openmmdl_protein_ligand_quickrun.py not found next to this script."
  exit 1
fi

export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"

if ! micromamba env list | awk '{print $1}' | grep -qx "openmmdl"; then
  echo "ERROR: environment 'openmmdl' not found. Run install script first."
  exit 1
fi
micromamba activate openmmdl

for f in "5wyz-moe-processed_openMMDL.pdb" "5VF.sdf"; do
  if [[ ! -f "${TUTORIAL_DIR}/${f}" ]]; then
    echo "ERROR: missing tutorial file ${TUTORIAL_DIR}/${f}"
    exit 1
  fi
done

cp "${SCRIPT_DIR}/openmmdl_protein_ligand_quickrun.py" "${OPENMMDL_DIR}/openmmdl_protein_ligand_quickrun.py"

cd "${OPENMMDL_DIR}"
rm -rf "${RUN_DIR}"
openmmdl simulation \
  -f "${RUN_DIR}" \
  -s "${OPENMMDL_DIR}/openmmdl_protein_ligand_quickrun.py" \
  -t "${TUTORIAL_DIR}/5wyz-moe-processed_openMMDL.pdb" \
  -l "${TUTORIAL_DIR}/5VF.sdf"

echo "Simulation finished. Output:"
ls -la "${RUN_DIR}"
