#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Legacy Backward-Compatibility Forwarder
# ==============================================================================
# This script forwards all arguments to the unified 'scripts/colab.sh'
# to ensure existing notebooks, bookmarks, and URLs continue to work seamlessly.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/colab.sh" "$@"
