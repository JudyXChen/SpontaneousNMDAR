#!/usr/bin/env bash
# Missing calcium_total runs for the Nov 2 slides, most important first (26 runs):
#   slide 24 EPSP + bAP:  10 NMDAR/cluster at 50/100/200 nm, seeds 1 and 5
#   slide 23 one cluster: 10 NMDAR, seeds 1-5
#   slide 24 glutamate only: 10 NMDAR/cluster at 50/100/200 nm, seeds 1-3, then 4-5
# Output: results/calcium_total/nov2/{evoked,spont}/CASE/seed_NNNNN
# Finished runs are skipped, so the same command resumes after an interruption.
# Usage: bash run/run_nov2_missing.sh [--smoke] [JOBS, default 5]
# Elsewhere than this Mac, set MCELL_PATH and PYTHON first.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
source scripts/runtime.sh
[[ -f "${MCELL_PATH:-}/lib/mcell.so" ]] || { echo 'Set MCELL_PATH to the folder containing lib/mcell.so' >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

if [[ "${1:-}" == --one ]]; then
    protocol="$2"; case_name="$3"; seed="$4"
    work="$NOV2_OUT/$protocol/$case_name/seed_$(printf %05d "$seed")"
    [[ -e "$work/DONE" ]] && { echo "SKIP $protocol $case_name seed $seed"; exit 0; }
    [[ -e "$work" ]] && mv "$work" "$work.partial_$(date +%s)"
    mkdir -p "$(dirname "$work")"
    echo "START $protocol $case_name seed $seed $(date '+%F %T')"
    status=0
    if [[ "$protocol" == evoked ]]; then
        env MCELL_ITERATIONS="$NOV2_ITERATIONS" MCELL_TIME_STEP=1e-7 MCELL_GLU_SITES=both \
            "$PYTHON" -B model/shared/run_legacy.py "model/calcium_total/evoked/$case_name" "$work" "$seed" >/dev/null || status=$?
    else
        "$PYTHON" -B "model/calcium_total/spont/$case_name/model.py" --seeds "$seed" --jobs 1 --output "$work" \
            ${NOV2_SMOKE:+--smoke} >/dev/null || status=$?
    fi
    if [[ $status -eq 0 ]]; then
        touch "$work/DONE"
        echo "DONE  $protocol $case_name seed $seed $(date '+%F %T')"
    else
        echo "FAILED $protocol $case_name seed $seed (status $status); see $work" >&2
    fi
    exit 0
fi

smoke=''
if [[ "${1:-}" == --smoke ]]; then smoke=1; shift; fi
jobs="${1:-5}"
export NOV2_OUT="${NOV2_OUT:-$repo/results/calcium_total/nov2}"
export NOV2_ITERATIONS=350000
if [[ -n "$smoke" ]]; then export NOV2_ITERATIONS=1000 NOV2_SMOKE=1; fi
mkdir -p "$NOV2_OUT"

queue() {
    for seed in 1 5; do for d in 50 100 200; do echo "evoked 2_cluster_${d}nm $seed"; done; done
    for seed in 1 2 3 4 5; do echo "evoked 1_cluster_10_NMDAR $seed"; done
    for seed in 1 2 3 4 5; do for d in 50 100 200; do echo "spont 2_cluster_${d}nm $seed"; done; done
}
[[ -x /usr/bin/caffeinate ]] && keep_awake=(/usr/bin/caffeinate -i) || keep_awake=()
queue | "${keep_awake[@]+"${keep_awake[@]}"}" xargs -P "$jobs" -L 1 bash "$repo/run/run_nov2_missing.sh" --one 2>&1 | tee -a "$NOV2_OUT/progress.log"
echo "Finished: $(find "$NOV2_OUT" -name DONE | wc -l | tr -d ' ') of 26 runs. Failures, if any: grep FAILED $NOV2_OUT/progress.log"
