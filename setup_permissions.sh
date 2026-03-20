#!/usr/bin/env bash
set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${REPO_ROOT}/scripts"

if [[ ! -d "${SCRIPTS_DIR}" ]]; then
  echo "[ERROR] scripts directory not found: ${SCRIPTS_DIR}" >&2
  exit 1
fi

find "${SCRIPTS_DIR}" -type f \( -name "*.py" -o -name "*.sh" \) -exec chmod +x {} +

echo "[INFO] Executable permissions have been applied to .py and .sh files under:"
echo "       ${SCRIPTS_DIR}"
echo
echo "[INFO] To use the top-level pipeline directly in the current shell, run:"
echo "       export PATH=\"${SCRIPTS_DIR}:\$PATH\""
echo
echo "[INFO] To make this change persistent, add the following line to ~/.bashrc:"
echo "       export PATH=\"${SCRIPTS_DIR}:\$PATH\""