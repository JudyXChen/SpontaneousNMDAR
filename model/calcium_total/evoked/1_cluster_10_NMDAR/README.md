> Current layout: recorded data is in `results/calcium_total/evoked/conditions/1_cluster_10_NMDAR/`.
> For current run commands and shared settings, see [the shared model guide](../../../../doc/shared_model.md).
> The notes below describe the historical condition.

# Calibrated 1 mM NMDAR block test

This folder was copied from `archive/NMDAR_old/model/evoked/1_cluster_10_NMDAR/`, then
only `k_B(V)` and `k_U(V)` were calibrated. The original model leaf and its
seed-1 outputs remain in the old-kinetics branch.

## Spatial run

- Condition: one cluster, 10 NMDAR, 500 released glutamate molecules
- Seed: 1
- Timestep: `1e-7` s
- Duration: 350,000 iterations (35 ms)
- Magnesium representation: the model's implicit 1 mM rate-table reference
- Status: `FINISHED`

| Metric | Original rates | Calibrated rates |
| --- | ---: | ---: |
| Peak sampled open NMDAR | 1 | 1 |
| Samples with at least one open NMDAR | 32 | 385 |
| Sampled open runs | 7 | 29 |
| Integrated sampled open occupancy (receptor ms) | 0.320 | 3.850 |
| Cumulative modeled NMDAR calcium-entry events | 202 | 1,323 |

Open occupancy was sampled every 10 us. Consecutive positive samples are not
microscopic dwell-time measurements because faster block/unblock transitions
can occur between samples. The 12.03-fold occupancy increase and 6.55-fold
calcium-entry increase describe one matched stochastic seed, not an ensemble
mean. The run also reported missed IP3R reactions (4.536% and 11.876%) and a
small SERCA miss rate (0.003%); no missed NMDAR reaction was reported.

The comparison data, metrics, and title-free plots are in
`archive/NMDAR_old/comparison/evoked/compare_1_cluster_10_NMDAR_seed1/`.
Raw MCell outputs are in `results/calcium_total/evoked/conditions/1_cluster_10_NMDAR`, under `react_data/`,
`viz_data/`, `logs/` and `reports/`.

## Isolated kinetic test

`archive/NMDAR_old/comparison/kinetics/100uM_Glu_1ms/`
contains an exact ten-state calculation at -67 mV and 1 mM
Mg for a controlled 100 uM, 1 ms glutamate pulse. Calibration increased peak
expected open occupancy for 10 receptors from 0.0228 to 0.0759 and integrated
open occupancy from 0.516 to 1.786 receptor ms. This isolates the kinetic-rate
effect but does not reproduce the spatial cleft glutamate waveform.

Neither test directly reproduces Espinosa and Kavalali's experiment, which
used voltage clamp, 1.25 or 2 mM Mg, approximately 100-event averages, and
reported the AP-5-sensitive fraction of total AMPAR-plus-NMDAR mEPSC charge.
