# Results

Simulation output is not stored in Git because a single run can reach about 1 GB.
Running the launchers fills this folder locally:

```text
results/
  calcium_total/
    evoked/runs/CASE/seed_NNNNN/  every evoked (EPSP + bAP) run: react_data, viz_data, logs, reports
    evoked/runs/CASE/plots/       figures for that one condition
    evoked/compare/NAME/          figures comparing conditions
    evoked/batches/               timestamped output of run_evoked_distance.sh; move finished seeds into runs/
    spont/runs, spont/compare     the same, for resting voltage (glutamate only)
    conductance/                  NMDAR conductance calibration
  calcium_track/
    evoked/CASE/seed_NNNNN/       counts, frames, voxels, validation, audit
    evoked/CASE/compare_seeds/    across-seed summary
    evoked/compare_*/             across-distance comparison (older runs without ER)
    clamped/TEST/CASE/seed_NNNNN/ clamped runs with ER
```

CASE uses the model folder name, `<clusters>x<NMDAR per cluster>[_<distance>nm]_<glutamate per site>g`,
for example `1x10_500g` or `2x10_100nm_250g`; `_Aonly` marks glutamate released at cluster A only.

Each run's `input_reference.json` points to the frozen inputs that produced it.
