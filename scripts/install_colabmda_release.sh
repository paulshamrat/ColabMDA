#!/usr/bin/env bash
set -euo pipefail

# Install ColabMDA from GitHub Release assets without cloning the full repository.
#
# Usage:
#   bash scripts/install_colabmda_release.sh
#   bash scripts/install_colabmda_release.sh <tag> <install_dir>
#
# Examples:
#   bash scripts/install_colabmda_release.sh latest /content/colabmda
#   bash scripts/install_colabmda_release.sh v0.1.0 /content/colabmda

REPO="${COLABMDA_REPO:-paulshamrat/ColabMDA}"
TAG="${1:-latest}"
INSTALL_DIR="${2:-/content/colabmda}"

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

WHEEL_URL="$(
python3 - "${REPO}" "${TAG}" <<'PY'
import json
import sys
import urllib.request

repo = sys.argv[1]
tag = sys.argv[2]
try:
    if tag == "latest":
        api = f"https://api.github.com/repos/{repo}/releases/latest"
    else:
        api = f"https://api.github.com/repos/{repo}/releases/tags/{tag}"

    with urllib.request.urlopen(api, timeout=15) as r:
        data = json.load(r)

    assets = data.get("assets", [])
    wheels = [a["browser_download_url"] for a in assets if a.get("name", "").endswith(".whl")]
    print(wheels[0] if wheels else "")
except Exception:
    print("")
PY
)"

echo "[INFO] Installing into current Python environment"
python3 -m pip install --upgrade pip

if [[ -n "${WHEEL_URL}" ]]; then
  echo "[INFO] Downloading release wheel: ${WHEEL_URL}"
  if curl -fsSL "${WHEEL_URL}" -o colabmda.whl; then
    python3 -m pip install --upgrade ./colabmda.whl
    echo "[OK] Installed from GitHub Release wheel"
  else
    echo "[WARN] Release wheel download failed."
    WHEEL_URL=""
  fi
fi

if ! command -v colabmda >/dev/null 2>&1; then
  echo "[INFO] Fallback: trying PyPI package"
  if ! python3 -m pip install --upgrade colabmda; then
    echo "[WARN] PyPI install failed."
  fi
fi

if ! command -v colabmda >/dev/null 2>&1; then
  if [[ "${TAG}" == "latest" ]]; then
    REF="main"
  else
    REF="${TAG}"
  fi
  echo "[INFO] Fallback: installing from GitHub source (${REF})"
  python3 -m pip install --upgrade "git+https://github.com/${REPO}.git@${REF}"
fi

echo "[OK] ColabMDA installed from GitHub Release in ${INSTALL_DIR}"
echo "[OK] Try: colabmda --help"
