#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: $PYTHON_BIN was not found. Install Python 3.11 in WSL first." >&2
  exit 1
fi

"$PYTHON_BIN" -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e ".[dev]"
python -m pytest -q

echo
echo "Environment ready. In VS Code select: .venv/bin/python"
echo "Next command: cml-repairnet doctor"
