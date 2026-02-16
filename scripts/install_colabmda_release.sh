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
if tag == "latest":
    api = f"https://api.github.com/repos/{repo}/releases/latest"
else:
    api = f"https://api.github.com/repos/{repo}/releases/tags/{tag}"

with urllib.request.urlopen(api) as r:
    data = json.load(r)

assets = data.get("assets", [])
wheels = [a["browser_download_url"] for a in assets if a.get("name", "").endswith(".whl")]
if not wheels:
    raise SystemExit("ERROR: no wheel asset found in release.")
print(wheels[0])
PY
)"

echo "[INFO] Downloading: ${WHEEL_URL}"
curl -fL "${WHEEL_URL}" -o colabmda.whl

echo "[INFO] Installing into current Python environment"
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade ./colabmda.whl

echo "[OK] ColabMDA installed from GitHub Release in ${INSTALL_DIR}"
echo "[OK] Try: colabmda --help"
