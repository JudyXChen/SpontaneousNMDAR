#!/usr/bin/env bash
# Run the model for a list of seeds, sequentially, teeing each to logs/seed_N.log.
#
# The upstream SpineApparatus run_seeds.sh was an empty file, so this is new.
#
# Usage:
#   ./run_seeds.sh [SEED ...]
#
# Seeds default to "1". Iteration count and timestep are taken from the
# environment (customization.py reads MCELL_ITERATIONS / MCELL_TIME_STEP /
# MCELL_SEED), so parameters.py stays at its 1 s production defaults.
#
#   MCELL_ITERATIONS=100000 ./run_seeds.sh 1        # 10 ms smoketest at 1e-7
#   MCELL_ITERATIONS=6000000 ./run_seeds.sh 1       # 600 ms formal run
#   MCELL_ITERATIONS=6000000 ./run_seeds.sh 1 2 3   # three seeds, one after another
#
# NOTE: do not raise MCELL_TIME_STEP above 1e-7. The IP3-binding reactions already
# exceed reaction probability 1 at 1e-7 (~0.5% of firings missed) and the 24-state
# IP3R triples them. A coarser step suppresses IP3R opening as a numerical artifact
# -- which is indistinguishable from the result this model exists to measure.
# Keep MCELL_TIME_STEP at 1e-7 s when using these models.
set -euo pipefail

cd "$(dirname "$0")"

if [[ -z "${MCELL_PATH:-}" ]]; then
    echo "Error: MCELL_PATH is not set. It must point at the MCell4 install, e.g." >&2
    echo "  export MCELL_PATH=/Applications/Blender-4.4.3-CellBlender/blender.app/Contents/Resources/4.4/extensions/system/cellblender/extensions/mcell" >&2
    exit 1
fi

SEEDS=("$@")
if [[ ${#SEEDS[@]} -eq 0 ]]; then
    SEEDS=(1)
fi

PYTHON="${PYTHON:-python}"
mkdir -p logs

echo "MCELL_PATH   = $MCELL_PATH"
echo "ITERATIONS   = ${MCELL_ITERATIONS:-<parameters.py default: 10000000>}"
echo "TIME_STEP    = ${MCELL_TIME_STEP:-<parameters.py default: 1e-7>}"
echo "seeds        = ${SEEDS[*]}"
echo

for seed in "${SEEDS[@]}"; do
    log="logs/seed_${seed}.log"
    echo "=== seed $seed -> $log ==="
    # model.py must run with the model directory as CWD: the react_data/, logs/ and
    # viz_data/ output paths in observables.py are CWD-relative.
    "$PYTHON" model.py -seed "$seed" 2>&1 | tee "$log"
done

echo
echo "Done. Logs in logs/, counts in react_data/seed_*/."
