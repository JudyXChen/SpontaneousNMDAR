#!/usr/bin/env bash
# Sourced by the three suite launchers; explicit server overrides remain possible.
set -euo pipefail
resources=/Applications/Blender-4.4.3-CellBlender/Blender.app/Contents/Resources/4.4
if [[ -z "${MCELL_PATH:-}" && -f "$resources/extensions/system/cellblender/extensions/mcell/lib/mcell.so" ]]; then
    export MCELL_PATH="$resources/extensions/system/cellblender/extensions/mcell"
fi
if [[ -z "${PYTHON:-}" ]]; then
    selected_mcell="${MCELL_PATH:-}"
    if [[ -n "$selected_mcell" && "$selected_mcell/lib/mcell.so" -ef "$resources/extensions/system/cellblender/extensions/mcell/lib/mcell.so" ]]; then
        export PYTHON="$resources/python/bin/python3.11"
    else
        export PYTHON=python3
    fi
fi
