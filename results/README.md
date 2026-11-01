# Results

Simulation output is not stored in Git because a single run can reach about 1 GB.
Running the launchers fills this folder locally:

```text
results/
  calcium_total/
    evoked/conditions/CASE/       one condition, seed 1: data and figures
    evoked/distance_seeds2_4/CASE/ seeds 2-4 for the distance series
    evoked/compare_*/             comparison figures across conditions
    spont/conditions/CASE/        resting-voltage conditions
    spont/compare_*/              comparison figures
    conductance/                  NMDAR conductance calibration
  calcium_track/
    evoked/CASE/seed_NNNNN/       counts, frames, voxels, validation, audit
    evoked/CASE/compare_seeds/    across-seed summary
    evoked/compare_*/             across-distance comparison (older runs without ER)
    clamped/TEST/CASE/seed_NNNNN/ clamped runs with ER
```

Each run's `input_reference.json` points to the frozen inputs that produced it.
