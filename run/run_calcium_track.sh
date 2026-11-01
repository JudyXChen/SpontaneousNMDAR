#!/usr/bin/env bash
# Scope the sleep assertion to this job, including its postprocessing.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
source_dir="$repo/model/calcium_track"
cd "$repo"
resources=/Applications/Blender-4.4.3-CellBlender/Blender.app/Contents/Resources/4.4
mcell_path="$resources/extensions/system/cellblender/extensions/mcell"
python_bin="$resources/python/bin/python3.11"
if [[ "${1:-}" == --batch && $# -eq 3 ]]; then
    plan="$2"
    session="$3"
    [[ -f "$plan" ]]
    if tmux has-session -t "$session" 2>/dev/null; then
        echo "Session already exists: $session" >&2
        exit 1
    fi
    batch_script="$source_dir/batch.py"
    printf -v launch_command 'cd %q && exec /usr/bin/caffeinate -i /usr/bin/env MCELL_PATH=%q MPLCONFIGDIR=%q XDG_CACHE_HOME=%q %q -B %q %q' \
        "$repo" "$mcell_path" /private/tmp/calcium-mpl /private/tmp/calcium-cache "$python_bin" "$batch_script" "$plan"
elif [[ $# -eq 3 || $# -eq 4 ]]; then
    inputs="$1"
    output="$2"
    session="$3"
    seed="${4:-1}"
    [[ -f "$inputs/model.py" && ! -e "$output" ]]
    printf -v launch_command 'cd %q && exec /usr/bin/caffeinate -i /usr/bin/env MCELL_PATH=%q %q -B %q %q --inputs %q --seed %q' \
        "$repo" "$mcell_path" "$python_bin" "$source_dir/run.py" "$output" "$inputs" "$seed"
else
    echo 'usage: run_calcium_track.sh INPUT_DIRECTORY OUTPUT_DIRECTORY SESSION [SEED]' >&2
    echo '   or: run_calcium_track.sh --batch PLAN SESSION' >&2
    exit 2
fi
tmux new-session -d -s "$session" "$launch_command"
echo "Started $session under job-scoped caffeinate; no system power settings changed."
