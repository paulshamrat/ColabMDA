#!/usr/bin/env bash
# Activate the installed ColabMDA environment before running these scripts.
set -euo pipefail
PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
: "${COLABMDA_RUN_ID:?Set a new run ID, or the original ID to resume}"
[[ "$COLABMDA_RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || { echo 'Invalid run ID' >&2; exit 2; }
export PROJECT_ROOT COLABMDA_RUN_ID
export COLABMDA_PLATFORM=CUDA COLABMDA_REQUIRE_GPU=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
SYSTEM_DIR="$PROJECT_ROOT/data/02_simulation/$COLABMDA_RUN_ID/apo/WT"
export SYSTEM_DIR
mkdir -p "$PROJECT_ROOT/logs/02_simulation/$COLABMDA_RUN_ID"

verify_run() {
    python -c 'import os; from colabmda.core.examples import record_run; record_run(os.environ["PROJECT_ROOT"], os.environ["COLABMDA_RUN_ID"])'
}
