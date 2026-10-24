#!/usr/bin/env bash
#SBATCH --job-name=cmda-pilot
#SBATCH --partition=work1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --signal=INT@120
set -euo pipefail

# A thin scheduler adapter. Scientific implementation stays in cmda.
# Arguments: PROJECT_ROOT RUN_ID STATE VARIANT STAGE [REPLICA] [TOTAL_NS]
[[ $# -ge 5 && $# -le 7 ]] || { echo "Expected root run state variant eq|md [r1] [0.1]" >&2; exit 2; }
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
[[ "$stage" = eq || "$stage" = md ]] || { echo "Stage must be eq or md" >&2; exit 2; }
: "${COLABMDA_CONDA_ROOT:?Export the installed writable conda root}"
source "$COLABMDA_CONDA_ROOT/etc/profile.d/conda.sh"
conda activate colabmda
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK:-4}"
cd "$PROJECT_ROOT/data/$state/simulation/$COLABMDA_RUN_ID/$variant"
test -s 01_prep/system.prmtop
test -s 01_prep/system.inpcrd
test -s 01_prep/PREP_DONE
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
if [[ "$stage" = eq ]]; then
  cd 02_eq
  srun cmda sim --state "$state" --variant "$variant" --stage eq --protocol staged --root "$PROJECT_ROOT"
else
  cd 03_md
  srun cmda sim --state "$state" --variant "$variant" --stage md --protocol staged \
    --replica "$replica" --total-ns "$total_ns" --checkpoint-ps 100 --traj-interval 10 --root "$PROJECT_ROOT"
fi
