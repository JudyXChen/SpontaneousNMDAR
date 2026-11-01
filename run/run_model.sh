#!/usr/bin/env bash
# Run an exported evoked model into a new result directory.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/scripts/runtime.sh"
exec "$PYTHON" -B "$ROOT/model/shared/run_legacy.py" "$@"
