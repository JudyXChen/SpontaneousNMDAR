#!/usr/bin/env bash
# EPSP + bAP run: 200 nm, 1000 glutamate per site, original membrane buffer, ER, seeds 1-5.
# To add the Bartol buffer, add cytosolic_buffer:2_cluster_200nm_cytosolic_buffer_evoked to the loop.
# Usage: bash run/run_evoked_buffer_test.sh [SEEDS_AT_ONCE, default 5]
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
source scripts/runtime.sh
workers="${1:-5}"
track=model/calcium_track
for case in original_buffer:2_cluster_200nm_evoked; do
    name="${case%%:*}"; condition="${case##*:}"
    inputs="$track/snapshots/evoked_buffer_test_$name"
    [[ -e "$inputs" ]] || "$PYTHON" -B "$track/build.py" "$inputs" \
        --condition "$track/conditions/$condition.json" --glutamate 1000 1000
    bash run/run_calcium_track.sh --batch "$track/batches/evoked_buffer_test_$name.json" "evoked-$name" "$workers"
done
