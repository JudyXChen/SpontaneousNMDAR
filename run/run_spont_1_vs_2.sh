#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/scripts/runtime.sh"
exec "$PYTHON" -B -u "$ROOT/model/shared/batch.py" --suite spont_1_vs_2 "$@"
