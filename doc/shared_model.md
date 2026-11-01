# Shared runtime and spontaneous conditions

`model/calcium_total/spont/CASE/` holds two small files. `instantiation.json` selects
PSD triangles, receptor counts, and glutamate counts and times. `parameters.json` sets
duration, waveform onset, timestep and VGCC populations. The family's `channels.json`
defines the channel drive. `model/shared/build.py` combines these with the shared
chemistry into a self-contained snapshot.

## Protocol

Everything is placed at t = 0, with NMDARs placed before other surface molecules.
Glutamate is released at 10 ms (500 at one cluster, or 250 at each of two), and the
run stops at 45 ms with a 1e-7 s timestep. Spontaneous conditions hold channel rates
at −67 mV and 1 mM Mg²⁺.

## Run

```sh
bash run/run_spont_1_vs_2.sh
bash run/run_spont_distance.sh --seeds 1 2 3 4 --jobs 3
bash run/run_spont_distance.sh --list
bash run/run_spont_distance.sh --prepare-only --seeds 1
bash run/run_spont_distance.sh --smoke --seeds 1
```

`--smoke` shortens the timeline to 0.1 ms + 0.1 ms and only checks startup and output
arithmetic. Existing output directories are never overwritten.

## Outputs

Each call creates a new results directory containing `manifest.json` (protocol, seeds,
runtime, input hashes) and `status.tsv`/`status.json`.

| File | Contents |
|---|---|
| `sNMDAR_World.dat` | Open receptors |
| `sNMDAR_total_World.dat` | Receptors across all 10 states (conserved) |
| `sNMDAR_Ca_flux_World.dat` | Cumulative modeled calcium-entry firings |
| `NMDAR_*`, `VGCC_*` | Aggregates across populations |

Counts are sampled every 10 µs. `validate_run.py` checks population sums, conservation,
monotonic entry, and release sites and times. Average seed-level values across
independent runs, not time samples.
