#!/usr/bin/env bash
# Start the local batch in a persistent tmux session, with idle sleep inhibited.
# Usage: bash start_evoked_local.sh [options for run_evoked_distance.sh]
# Example short background test: bash start_evoked_local.sh --smoke
set -euo pipefail
TASK_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

# This entry point runs inside tmux; the public invocation below starts it.
if [[ "${1:-}" == --worker ]]; then
    shift
    run_output="$1"
    shift
    cd "$TASK_ROOT"
    # A tmux login shell can inherit a different PYTHON or MCELL_PATH.
    # This local-only launcher pins the runtime used by the successful smoke tests.
    local_resources=/Applications/Blender-4.4.3-CellBlender/Blender.app/Contents/Resources/4.4
    export MCELL_PATH="$local_resources/extensions/system/cellblender/extensions/mcell"
    export PYTHON="$local_resources/python/bin/python3.11"
    set +e
    /usr/bin/caffeinate -i bash "$TASK_ROOT/run/run_evoked_distance.sh" "$@" --output "$run_output" 2>&1 | tee "$run_output.launcher.log"
    batch_status=${PIPESTATUS[0]}
    printf '\nBatch exit status: %s\n' "$batch_status" | tee -a "$run_output.launcher.log"
    exit "$batch_status"
fi

command -v tmux >/dev/null || { echo 'tmux is required.' >&2; exit 1; }
[[ -x /usr/bin/caffeinate ]] || { echo 'This launcher requires macOS caffeinate.' >&2; exit 1; }
session=evoked-distance
if tmux has-session -t "=$session" 2>/dev/null; then
    echo "Session $session already exists; no new batch started."
    echo "View it with: tmux attach -t $session"
    exit 1
fi
run_parent="$TASK_ROOT/results/calcium_total/evoked/distance_seeds2_4"
mkdir -p "$run_parent"
run_output="$run_parent/local_$(date -u +%Y%m%dT%H%M%SZ)_$$"
printf -v worker_command '%q ' bash "$TASK_ROOT/run/start_evoked_local.sh" --worker "$run_output" "$@"
tmux new-session -d -s "$session" "$worker_command"
printf '%s\n' "$run_output" > "$run_parent/local_latest.txt"
echo "Started tmux session: $session"
echo "Output directory: $run_output"
echo "Launcher log: $run_output.launcher.log"
echo "View live: tmux attach -t $session"
echo 'Detach: Ctrl+b, release, then d. Keep the laptop plugged in with the lid open.'
