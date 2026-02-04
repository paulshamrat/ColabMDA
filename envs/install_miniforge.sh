#!/usr/bin/env bash
set -euo pipefail

PREFIX=${MINIFORGE_PREFIX:-$HOME/miniforge3}

wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
if [[ -d "$PREFIX" ]]; then
  bash /tmp/miniforge.sh -b -u -p "$PREFIX"
else
  bash /tmp/miniforge.sh -b -p "$PREFIX"
fi
export PATH="$PREFIX/bin:$PATH"
# shellcheck source=/dev/null
source "$PREFIX/etc/profile.d/conda.sh"

conda --version
