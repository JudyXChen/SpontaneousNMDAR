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
| Glutamate | Per site; default 250, released at t = 0. `build.py --glutamate 1000 0` releases 1000 at A only |
| Voltage | Clamped: −67 mV, 1 mM Mg²⁺, so all channel rates are constant |
| Duration | 35 ms; pumps and buffers start unbound, with no equilibration |
| Timestep | 0.1 µs; counts and particle frames every 10 µs |
| Geometry | Spine and dendrite with the smooth, unlaminated ER (`SA` mesh), plus SERCA, 60 RyR and IP3R |
| Immobile buffer | `membrane`: 4791/µm² on the spine membrane, kon 10⁶ M⁻¹s⁻¹, koff 2 s⁻¹. `cytosol`: 78.7 µM in the spine cytoplasm, kon 2.47×10⁸ M⁻¹s⁻¹, koff 524 s⁻¹ (Bartol et al. 2015) |

## Chemistry

Labels persist through mobile and fixed buffers, pumps and PLC. Unoccupied partners
are shared, so all calcium species compete for the same capacity. Extrusion removes
the label; pump leaks produce unlabeled `Ca2`.

Labels stop at the ER. Labeled calcium binds SERCA, RyR and IP3R and enters through
open RyRs at the same rates as `Ca2`, and that step ends its label. Calcium returning
from the ER or from these proteins is unlabeled `Ca2`.

## Outputs

| File | Contents |
|---|---|
| `react_data/.../NMDAR_A_Ca_flux_World.dat` | Cumulative entry through cluster A (B likewise) |
| `react_data/.../Ca_A_World.dat` | Free `Ca_A` |
| `react_data/.../Ca_A_inventory_World.dat` | Free plus bound `Ca_A` |
| `react_data/.../Ca_A_to_RyR_World.dat` | Cumulative `Ca_A` ending at RyR (also `_to_SERCA`, `_to_IP3R`) |
| `react_data/.../RyR_O_World.dat`, `Ca2_ER_SA.dat` | Open RyRs; ER calcium |
| `validation.json`, `audit.json` | Time grids, source balance, capacity and placement checks |
| `spatial/calcium_voxels.csv.gz` | Sparse counts by time, voxel, source and bound state |
| `blender_ca_only/` | Clean movie of free `Ca_A` and `Ca_B` |

Validation requires `entry = free + bound + removed + handed to the ER` for each label at every count
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
