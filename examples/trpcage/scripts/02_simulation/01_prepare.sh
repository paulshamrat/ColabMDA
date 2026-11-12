#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
check_stage "$PROJECT_ROOT/data/01_structure/trpcage" structure
verify_run
cmda prep --root "$PROJECT_ROOT" --state apo --variant WT --receptor "$PROJECT_ROOT/data/01_structure/trpcage/WT.pdb"
for name in model.pdb solvated.pdb system.prmtop system.inpcrd; do test -s "$SYSTEM_DIR/01_prep/$name"; done
complete_stage "$SYSTEM_DIR/01_prep" prep "Preparation outputs verified"
