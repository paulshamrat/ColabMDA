#!/usr/bin/env bash
#SBATCH --job-name=cmda-pilot
#SBATCH --partition=work1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --signal=INT@120
set -euo pipefail

# A thin scheduler adapter. Scientific implementation stays in cmda.
# Arguments: PROJECT_ROOT RUN_ID STATE VARIANT STAGE [REPLICA] [TOTAL_NS]
[[ $# -ge 5 && $# -le 7 ]] || { echo "Expected root run state variant eq|md|smoke [r1] [0.1]" >&2; exit 2; }
export PROJECT_ROOT="$1"
export COLABMDA_RUN_ID="$2"
state="$3"
variant="$4"
stage="$5"
replica="${6:-r1}"
total_ns="${7:-0.1}"
[[ "$PROJECT_ROOT" = /* && -d "$PROJECT_ROOT" ]] || { echo "Use an existing absolute scratch project" >&2; exit 2; }
for value in "$COLABMDA_RUN_ID" "$state" "$variant" "$replica"; do
  [[ "$value" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || { echo "Invalid run/state/variant/replica" >&2; exit 2; }
done
[[ "$state" = apo || "$state" = atpmg ]] || { echo "This adapter supports apo or atpmg" >&2; exit 2; }
[[ "$stage" = eq || "$stage" = md || "$stage" = smoke ]] || { echo "Stage must be eq, md or smoke" >&2; exit 2; }
: "${COLABMDA_CONDA_ROOT:?Export the installed conda root}"
source "$COLABMDA_CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "${COLABMDA_ENV:-colabmda}"
export COLABMDA_REQUIRE_GPU=1
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK:-4}"
cd "$PROJECT_ROOT/data/$state/simulation/$COLABMDA_RUN_ID/$variant"
test -s 01_prep/system.prmtop
test -s 01_prep/system.inpcrd
test -s 01_prep/PREP_DONE
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
if [[ "$stage" = smoke ]]; then
  srun cmda sim --state "$state" --variant "$variant" --stage em --protocol staged --root "$PROJECT_ROOT"
  srun cmda sim --state "$state" --variant "$variant" --stage nvt,npt --protocol quick \
    --software-test --steps 50 --root "$PROJECT_ROOT"
  srun cmda sim --state "$state" --variant "$variant" --stage md --protocol quick \
    --software-test --steps 50 --checkpoint-ps 0.1 --replica r1 --root "$PROJECT_ROOT"
  python - <<'PY_CHECK'
from pathlib import Path
import MDAnalysis as mda
import numpy as np
trajectory = Path("03_md/r1/prod_0to0.1ps.dcd")
universe = mda.Universe("01_prep/system.prmtop", str(trajectory))
assert len(universe.trajectory) == 5
for frame in universe.trajectory:
    assert np.isfinite(universe.atoms.positions).all()
assert Path("03_md/r1/prod.chk").stat().st_size > 0
Path("SMOKE_DONE").write_text("CUDA software smoke: five finite frames; not research equilibration\n")
print("GPU_SMOKE_DONE", len(universe.atoms), "atoms, five finite frames")
PY_CHECK
elif [[ "$stage" = eq ]]; then
  cd 02_eq
  srun cmda sim --state "$state" --variant "$variant" --stage eq --protocol staged --root "$PROJECT_ROOT"
else
  cd 03_md
  srun cmda sim --state "$state" --variant "$variant" --stage md --protocol staged \
    --replica "$replica" --total-ns "$total_ns" --checkpoint-ps 100 --traj-interval 10 --root "$PROJECT_ROOT"
fi
