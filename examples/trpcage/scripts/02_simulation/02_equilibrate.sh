#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
# Short software check, not a production equilibration recommendation.
EQUILIBRATION_PS=100
test -s "$SYSTEM_DIR/01_prep/PREP_DONE"
verify_run
cmda sim --root "$PROJECT_ROOT" --state apo --variant WT --stage eq --protocol quick --equil-time "$EQUILIBRATION_PS"
for name in system.xml npt.state.xml equilibration_qc.json; do test -s "$SYSTEM_DIR/02_eq/$name"; done
printf 'Quick protocol outputs verified; inspect physical QC separately\n' > "$SYSTEM_DIR/02_eq/EQ_DONE"
