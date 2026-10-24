#!/usr/bin/env bash
# A small public-CLI workflow; all runtime products belong to the supplied project.
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo 'Usage: bash run.sh PROJECT_ROOT RECEPTOR_PDB [prepare|equilibrate|simulate|analyze|all]' >&2
    exit 2
fi
: "${COLABMDA_RUN_ID:?Set a fresh run ID, or retain the original ID to resume}"
project_root="$1"
receptor_path="$2"
stage="${3:-all}"
state="${COLABMDA_STATE:-apo}"
variant="${COLABMDA_VARIANT:-WT}"
export COLABMDA_PLATFORM="${COLABMDA_PLATFORM:-CUDA}"
export COLABMDA_REQUIRE_GPU=1
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
case "$stage" in
    prepare|equilibrate|simulate|analyze|all) ;;
    *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
if [[ "$stage" == prepare || "$stage" == all ]]; then
    cmda prep --root "$project_root" --state "$state" --variant "$variant" --receptor "$receptor_path"
fi
if [[ "$stage" == equilibrate || "$stage" == all ]]; then
    cmda sim --root "$project_root" --state "$state" --variant "$variant" --stage eq --protocol quick --equil-time 100
fi
if [[ "$stage" == simulate || "$stage" == all ]]; then
    cmda sim --root "$project_root" --state "$state" --variant "$variant" --stage md --replica r1 --total-ns "${COLABMDA_EXAMPLE_NS:-1}" --checkpoint-ps 100 --traj-interval 10
fi
if [[ "$stage" == analyze || "$stage" == all ]]; then
    cmda merge --root "$project_root" --state "$state" --variant "$variant" --replicas r1
    cmda analyze --root "$project_root" --state "$state" --variant "$variant"
fi
