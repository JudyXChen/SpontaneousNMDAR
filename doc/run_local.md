Local evoked distance batch

The defaults are experiments 2 and 3, each at 50/100/200 nm, with seeds
2, 3 and 4. This is 18 simulations, with three running concurrently.
Each simulation lasts 35 ms (350000 iterations at 1e-7 s).
The six distance models place NMDARs before the other released molecules
to prevent their small cluster regions from filling before placement.
Counts, geometry, and release times are unchanged; seeded trajectories
differ from the earlier initialization order.
The installed CellBlender Python/MCell paths are selected automatically.

Start from the project directory:
  bash run/start_evoked_local.sh

This creates the detached tmux session evoked-distance and returns your
terminal immediately. You may close Terminal or disconnect Codex; tmux
keeps the batch running. Do not start a second copy while it is running.

Watch the live progress:
  tmux attach -t evoked-distance

Detach again without stopping anything: Ctrl+b, release, then d.
To deliberately stop the batch, attach and press Ctrl+c. The launcher
stops its active simulations and records the interrupted runs.

Check results without attaching (run these from the project directory):
  run_output=$(cat results/calcium_total/evoked/distance_seeds2_4/local_latest.txt)
  cat "$run_output/status.tsv"
  tail -20 "$run_output.launcher.log"

The tmux session closes when the batch finishes. The log reports its final
exit status: 0 means all simulations passed the output checks. Nonzero
means failure/interruption; inspect status.json and the affected run.log.

Outputs are isolated under results/calcium_total/evoked/distance_seeds2_4/local_TIMESTAMP_PID/.
Existing seed-1 outputs and model input files are preserved. Each run
contains an input reference, manifest, logs, counts and visualization
files. Exact input snapshots live under model/shared/snapshots/.
This launcher does not produce multiseed comparison plots.

The launcher uses caffeinate -i to prevent idle system sleep while it runs.
Keep the Mac plugged in with the lid open. The display may turn off.
This does not make simulations continue through shutdown or lid sleep.

Optional short background test (1000 iterations per simulation):
  bash run/start_evoked_local.sh --smoke

To reduce concurrent simulations further, while retaining all three seeds:
  bash run/start_evoked_local.sh --jobs 2

The separate doc/run_remote.md guide remains available for future SSH
use. Its five-seed example supplies the seed list and worker count explicitly.
