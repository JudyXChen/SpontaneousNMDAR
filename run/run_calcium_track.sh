#!/usr/bin/env bash
# Start a calcium_track run or batch in a detached tmux session.
# Elsewhere than this Mac, set MCELL_PATH and PYTHON first.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
source_dir="$repo/model/calcium_track"
cd "$repo"
source "$repo/scripts/runtime.sh"
[[ -n "${MCELL_PATH:-}" && -f "$MCELL_PATH/lib/mcell.so" ]] || { echo 'Set MCELL_PATH to the folder containing lib/mcell.so' >&2; exit 1; }
command -v tmux >/dev/null || { echo 'tmux is required.' >&2; exit 1; }
keep_awake=()
[[ -x /usr/bin/caffeinate ]] && keep_awake=(/usr/bin/caffeinate -i)
cache="${TMPDIR:-/tmp}/calcium-track-cache"
mkdir -p "$cache"
if [[ "${1:-}" == --batch && ( $# -eq 3 || $# -eq 4 ) ]]; then
    plan="$2"; session="$3"; workers="${4:-}"
    [[ -f "$plan" ]]
    if tmux has-session -t "=$session" 2>/dev/null; then
        echo "Session already exists: $session" >&2
        exit 1
    fi
    command=("$PYTHON" -B "$source_dir/batch.py" "$plan")
    [[ -n "$workers" ]] && command+=(--workers "$workers")
elif [[ $# -eq 3 || $# -eq 4 ]]; then
    inputs="$1"; output="$2"; session="$3"; seed="${4:-1}"
    [[ -f "$inputs/model.py" && ! -e "$output" ]]
    command=("$PYTHON" -B "$source_dir/run.py" "$output" --inputs "$inputs" --seed "$seed")
else
    echo 'usage: run_calcium_track.sh INPUT_DIRECTORY OUTPUT_DIRECTORY SESSION [SEED]' >&2
    echo '   or: run_calcium_track.sh --batch PLAN SESSION [WORKERS]' >&2
    exit 2
fi
printf -v launch 'cd %q && exec' "$repo"
printf -v rest ' %q' "${keep_awake[@]+"${keep_awake[@]}"}" /usr/bin/env "MCELL_PATH=$MCELL_PATH" "MPLCONFIGDIR=$cache" "XDG_CACHE_HOME=$cache" "${command[@]}"
tmux new-session -d -s "$session" "$launch$rest"
echo "Started tmux session $session"
