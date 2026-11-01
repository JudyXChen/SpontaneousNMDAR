#!/usr/bin/env bash
# Buffer test: 200 nm, 1000 glutamate per site, clamped voltage, ER, release at t=0, seeds 1-5.
# Original membrane buffer and Bartol cytosolic buffer, each in its own tmux session.
# Usage: bash run/run_buffer_test.sh [SEEDS_AT_ONCE_PER_BUFFER, default 5]
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
source scripts/runtime.sh
workers="${1:-5}"
track=model/calcium_track
for case in original_buffer:2_cluster_200nm cytosolic_buffer:2_cluster_200nm_cytosolic_buffer; do
    name="${case%%:*}"; condition="${case##*:}"
    inputs="$track/snapshots/buffer_test_$name"
    [[ -e "$inputs" ]] || "$PYTHON" -B "$track/build.py" "$inputs" \
        --condition "$track/conditions/$condition.json" --glutamate 1000 1000
    bash run/run_calcium_track.sh --batch "$track/batches/buffer_test_$name.json" "buffer-$name" "$workers"
done
