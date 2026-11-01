# Source-labeled calcium

Calcium entering through NMDAR cluster A is `Ca_A` (pink) and through cluster B is
`Ca_B` (green). VGCC and pump-leak calcium is unlabeled `Ca2`. All three have the
same diffusion constant and reaction rates. Labels identify the **entry cluster**,
not the glutamate release site: glutamate is shared and can reach both clusters.

## Protocol

| Setting | Value |
|---|---|
| NMDARs | 10 + 10, immobile |
| Separation | 101.5, 200.5 or 405.2 nm, mirrored around the PSD center |
| Glutamate | Per site; default 250, released at t = 10 ms. `build.py --glutamate 1000 0` releases 1000 at A only |
| Voltage | Clamped: −67 mV, 1 mM Mg²⁺, so all channel rates are constant |
| Duration | 45 ms: 10 ms baseline + 35 ms response |
| Timestep | 0.1 µs; counts every 10 µs |
| Geometry | Spine and attached dendrite; spine apparatus and ER chemistry removed |

**The 10 ms baseline does not fully equilibrate the calcium system.** With the
current rates, PMCA reaches about 5% occupancy at 10 ms against an expected resting
value near 26%, and the buffers have time constants of about 0.5–1 s.

## Chemistry

Labels persist through mobile and fixed buffers, pumps and PLC. Unoccupied partners
are shared, so all calcium species compete for the same capacity. Extrusion removes
the label; pump leaks produce unlabeled `Ca2`.

## Outputs

| File | Contents |
|---|---|
| `react_data/.../NMDAR_A_Ca_flux_World.dat` | Cumulative entry through cluster A (B likewise) |
| `react_data/.../Ca_A_World.dat` | Free `Ca_A` |
| `react_data/.../Ca_A_inventory_World.dat` | Free plus bound `Ca_A` |
| `validation.json`, `audit.json` | Time grids, source balance, capacity and placement checks |
| `spatial/calcium_voxels.csv.gz` | Sparse counts by time, voxel, source and bound state |
| `blender_ca_only/` | Clean movie of free `Ca_A` and `Ca_B` |

Validation requires `entry = free + bound + removed` for each label at every count
sample. Voxel `i` spans `[i*size, (i+1)*size)` µm on each axis. Voxels hold counts,
not concentrations, because boundary voxels are only partly cytosol.

## Plots

```sh
python3 -B scripts/plot_calcium_voxels.py NEW_OUTPUT RUN_1 RUN_2 ...
python3 -B model/calcium_track/analysis/plot.py RUN_DIRECTORY
python3 -B model/calcium_track/analysis/compare_spatial.py RUN_A RUN_B NEW_OUTPUT
```

`plot_calcium_voxels.py` averages seeds and shows where each label's free ions spend time,
split by time since their own cluster conducted. Use the CellBlender Python, or any
Python with NumPy and Matplotlib.
