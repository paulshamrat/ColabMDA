#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
test -s "$SYSTEM_DIR/03_md/r1/MD_DONE"
cmda merge --root "$PROJECT_ROOT" --state apo --variant WT --replicas r1
cmda analyze --root "$PROJECT_ROOT" --state apo --variant WT
for name in rmsd.csv rg.csv rmsf.csv analysis.json; do test -s "$SYSTEM_DIR/04_analysis/$name"; done
printf 'Numerical analysis and downstream figures generated\n' > "$SYSTEM_DIR/04_analysis/ANALYSIS_DONE"
