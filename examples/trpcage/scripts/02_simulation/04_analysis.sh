#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
check_stage "$SYSTEM_DIR/03_md/r1" md
cmda merge --root "$PROJECT_ROOT" --state apo --variant WT --replicas r1
cmda analyze --root "$PROJECT_ROOT" --state apo --variant WT
for name in rmsd.csv rg.csv rmsf.csv; do test -s "$SYSTEM_DIR/04_analysis/$name"; done
complete_stage "$SYSTEM_DIR/04_analysis" analysis "Numerical analysis and downstream figures generated"
