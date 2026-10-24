#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
test -s "$PROJECT_ROOT/data/01_structure/trpcage/STRUCTURE_DONE"
verify_run
cmda prep --root "$PROJECT_ROOT" --state apo --variant WT --receptor "$PROJECT_ROOT/data/01_structure/trpcage/WT.pdb"
for name in model.pdb solvated.pdb system.prmtop system.inpcrd prep.json; do test -s "$SYSTEM_DIR/01_prep/$name"; done
printf 'Preparation outputs verified\n' > "$SYSTEM_DIR/01_prep/PREP_DONE"
