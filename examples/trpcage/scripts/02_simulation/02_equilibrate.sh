#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
# Short software check, not a production equilibration recommendation.
EQUILIBRATION_PS=100
check_stage "$SYSTEM_DIR/01_prep" prep
verify_run
cmda sim --software-test --root "$PROJECT_ROOT" --state apo --variant WT --stage eq --protocol quick --equil-time "$EQUILIBRATION_PS"
for name in system.xml npt.state.xml eq.pdb; do test -s "$SYSTEM_DIR/02_eq/$name"; done
complete_stage "$SYSTEM_DIR/02_eq" eq "Quick protocol outputs verified; inspect physical QC separately"
