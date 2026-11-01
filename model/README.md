# Model files

Each MCell model uses the same small set of files:

| File | Contents |
|---|---|
| `model.py` | Assembles the model and runs it |
| `parameters.py` | Seed, timestep, iterations and BNGL parameters |
| `subsystem.py` | Species, surface classes and reaction rules |
| `instantiation.py` | Geometry and release sites: where everything starts |
| `observables.py` | What is counted and written to disk |
| `customization.py` | Run-time options |
| `shared.py` | Holds parameter overrides |
| `geometry/` | Spine, dendrite and presynaptic meshes |
| `model.bngl` | Parameters and reactions in BNGL |

## Where each piece lives

| Folder | Role |
|---|---|
| `shared/geometry/` | The one copy of the meshes |
| `shared/background.bngl`, `nmda_kinetics.json` | Non-channel chemistry; NMDAR gating |
| `shared/resting_rates.json`, `evoked_rates.json` | Channel rates at −67 mV; EPSP + bAP tables |
| `shared/runtime/` | The files above, driven by small JSON conditions |
| `shared/legacy/` | The same files as originally exported from CellBlender; used by `calcium_total/evoked/` |
| `shared/build.py`, `batch.py`, `validate_run.py`, `run_legacy.py` | Build, run and check runs |
| `calcium_total/evoked/CASE/` | One folder per condition: own `model.bngl` and `instantiation.py`, links for the rest |
| `calcium_total/spont/CASE/` | One folder per condition: `instantiation.json` and `parameters.json` |
| `calcium_track/` | See below |

## calcium_track

| Path | Contents |
|---|---|
| `conditions/` | One JSON per condition: cluster triangles, receptors, glutamate, timing |
| `chemistry/` | Shared chemistry and the rules that label calcium by source |
| `runtime/` | Model files; files that differ from `shared/runtime/` are real copies, the rest are links |
| `geometry/` | Links to the shared meshes, including the smooth ER |
| `build.py` | Combines these into a frozen, self-contained input folder |
| `run.py`, `batch.py`, `batches/` | Run one seed or a planned list of seeds, then analyze |
| `analysis/` | Validation, voxel counts, audits and plots |

Conditions select PSD mesh triangles by index, so realized cluster distances are
close to, but not exactly, the name. See
[cluster geometry](../doc/cluster_geometry.md).

Links in this folder are relative symbolic links: clone with Git rather than
downloading a ZIP on Windows.
