#!/usr/bin/env bash
#SBATCH --job-name=cmda-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail

# CPU-only downstream processing; numerical and plotting engines stay in cmda.
[[ $# -eq 4 ]] || { echo "Expected root run state variant" >&2; exit 2; }
export PROJECT_ROOT="$1"
export COLABMDA_RUN_ID="$2"
state="$3"
variant="$4"
for value in "$COLABMDA_RUN_ID" "$state" "$variant"; do
  [[ "$value" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || { echo "Invalid run/state/variant" >&2; exit 2; }
done
[[ "$PROJECT_ROOT" = /* && -d "$PROJECT_ROOT" ]] || exit 2
[[ "$state" = apo || "$state" = atpmg ]] || exit 2
: "${COLABMDA_CONDA_ROOT:?Export the installed conda root}"
source "$COLABMDA_CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "${COLABMDA_ENV:-colabmda}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
run="$PROJECT_ROOT/data/$state/simulation/$COLABMDA_RUN_ID/$variant"
analysis="$PROJECT_ROOT/data/$state/analysis/$COLABMDA_RUN_ID/$variant/r1"
rm -f "$analysis/ANALYSIS_DONE"
cd "$run/03_md"
cmda merge --state "$state" --variant "$variant" --replicas r1 --root "$PROJECT_ROOT"
mkdir -p "$analysis"
cd "$analysis"
cmda analyze --state "$state" --variant "$variant" --replicas r1 --root "$PROJECT_ROOT"
for output in analysis.json rmsd.csv rg.csv rmsf.csv figures/rmsd_vs_time.png figures/rg_vs_time.png figures/rmsf_per_residue.png; do
  test -s "$output"
done
printf 'Verified r1 tables and figures generated\n' > ANALYSIS_DONE
