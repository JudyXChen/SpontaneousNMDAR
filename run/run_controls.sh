#!/usr/bin/env bash
# calcium_total evoked control experiments, seeds 1-5 (60 runs), most important first:
#   part 1, cluster split:   1x10 NMDAR / 250 glu;  1x20 NMDAR / 500 glu (area-matched patch)
#   part 2, glutamate dose:  2x10 NMDAR at 50/100/200 nm, 500 and 1000 glu per site
#   single-cluster dose:     1x10 NMDAR / 1000 glu
#   crosstalk:               2x10 NMDAR at 50/100/200 nm, 250 glu at cluster A only
# Baselines already run elsewhere (1x10/500; 2x10/250+250) are not repeated.
# Output: results/calcium_total/controls/LABEL/seed_NNNNN, where LABEL is the model
# folder name, plus _gluA for the crosstalk runs.
# Finished runs are skipped, so the same command resumes after an interruption.
# Usage: bash run/run_controls.sh [--smoke] [JOBS, default 8]
# Elsewhere than this Mac, set MCELL_PATH and PYTHON first.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
source scripts/runtime.sh
[[ -f "${MCELL_PATH:-}/lib/mcell.so" ]] || { echo 'Set MCELL_PATH to the folder containing lib/mcell.so' >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

if [[ "${1:-}" == --one ]]; then
    case_name="$2"; sites="$3"; seed="$4"
    label="$case_name"; [[ "$sites" == 1 ]] && label="${case_name}_gluA"
    work="$CONTROLS_OUT/$label/seed_$(printf %05d "$seed")"
    [[ -e "$work/DONE" ]] && { echo "SKIP  $label seed $seed"; exit 0; }
    [[ -e "$work" ]] && mv "$work" "$work.partial_$(date +%s)"
    mkdir -p "$(dirname "$work")"
    echo "START $label seed $seed $(date '+%F %T')"
    status=0
    env MCELL_ITERATIONS="$CONTROLS_ITERATIONS" MCELL_TIME_STEP=1e-7 MCELL_GLU_SITES="$sites" \
        "$PYTHON" -B model/shared/run_legacy.py "model/calcium_total/evoked/$case_name" "$work" "$seed" >/dev/null || status=$?
    if [[ $status -eq 0 ]]; then
        # Every receptor must be placed, and every trace must cover the full time grid.
        expected=$(awk '$1 == "n_NMDAR" {print $2}' "model/calcium_total/evoked/$case_name/model.bngl")
        placed=$(grep -o 'Released [0-9]* NMDAR_0' "$work/logs/seed_$seed.log" | awk '{s += $2} END {print s + 0}')
        [[ "$placed" == "$expected" ]] || { echo "placed $placed of $expected NMDAR" >&2; status=97; }
        rows=$((CONTROLS_ITERATIONS / 100 + 1))
        for name in NMDAR_World NMDAR_Ca_flux_World Ca2_Mushroom_transparent RyR_O_World; do
            n=$(wc -l < "$work/react_data/seed_$(printf %05d "$seed")/$name.dat" 2>/dev/null || echo 0)
            [[ $n -eq $rows ]] || { echo "$name.dat has $n of $rows rows" >&2; status=98; }
        done
    fi
    if [[ $status -eq 0 ]]; then
        touch "$work/DONE"
        echo "DONE  $label seed $seed $(date '+%F %T')"
    else
        echo "FAILED $label seed $seed (status $status); see $work" >&2
    fi
    exit 0
fi

smoke=''
if [[ "${1:-}" == --smoke ]]; then smoke=1; shift; fi
jobs="${1:-8}"
export CONTROLS_OUT="${CONTROLS_OUT:-$repo/results/calcium_total/controls}"
export CONTROLS_ITERATIONS=350000
if [[ -n "$smoke" ]]; then
    export CONTROLS_ITERATIONS=1000 CONTROLS_OUT="${CONTROLS_OUT}_smoke"
fi
mkdir -p "$CONTROLS_OUT"

queue() {
    for c in 1_cluster_10_NMDAR_250_glu 1_cluster_20_NMDAR_500_glu; do
        for seed in 1 2 3 4 5; do echo "$c both $seed"; done; done
    for g in 500 1000; do for d in 50 100 200; do
        for seed in 1 2 3 4 5; do echo "2_cluster_${d}nm_${g}_glu both $seed"; done; done; done
    for seed in 1 2 3 4 5; do echo "1_cluster_10_NMDAR_1000_glu both $seed"; done
    for d in 50 100 200; do for seed in 1 2 3 4 5; do echo "2_cluster_${d}nm 1 $seed"; done; done
}
total=$(queue | wc -l | tr -d ' ')
[[ -x /usr/bin/caffeinate ]] && keep_awake=(/usr/bin/caffeinate -i) || keep_awake=()
queue | "${keep_awake[@]+"${keep_awake[@]}"}" xargs -P "$jobs" -L 1 bash "$repo/run/run_controls.sh" --one 2>&1 | tee -a "$CONTROLS_OUT/progress.log"
echo "Finished: $(find "$CONTROLS_OUT" -name DONE | wc -l | tr -d ' ') of $total runs. Failures, if any: grep FAILED $CONTROLS_OUT/progress.log"
