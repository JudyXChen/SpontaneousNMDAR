# NMDAR clusters in a dendritic spine: MCell models

Stochastic MCell4 simulations of glutamate release, NMDAR activation and calcium
entry in a reconstructed mushroom spine, testing how far apart two NMDAR clusters
must sit to signal independently.

```text
model/
  shared/          common meshes, chemistry, rate tables and runtime
  calcium_total/   one pooled cytosolic calcium species (Ca2)
  calcium_track/   calcium labeled by entry cluster (Ca_A, Ca_B)
run/               simulation launchers
scripts/           figure scripts
doc/               model and method notes
results/           simulation output (see results/README.md)
```

## Requirements

MCell4 with its compatible Python (on Windows, see [doc/run_windows.md](doc/run_windows.md)), as shipped with
[CellBlender](https://mcell.org). Launchers find a local Blender 4.4.3 CellBlender
install automatically; elsewhere set `MCELL_PATH` and `PYTHON`.

## Run

From the repository root:

```sh
# Source-labeled calcium: build frozen inputs, then run one seed
python3 -B model/calcium_track/build.py model/calcium_track/snapshots/my_run \
  --condition model/calcium_track/conditions/2_cluster_200nm.json
bash run/run_calcium_track.sh model/calcium_track/snapshots/my_run \
  results/calcium_track/my_run my-session 1

# Immobile-buffer comparison: two batches of five seeds in tmux
bash run/run_buffer_test.sh

# The same comparison with EPSP + bAP
bash run/run_evoked_buffer_test.sh

# Pooled calcium at resting voltage
bash run/run_spont_distance.sh --list
bash run/run_spont_distance.sh --seeds 1 2 3 4

# Pooled calcium with the prescribed EPSP + bAP waveform
bash run/run_evoked_distance.sh --help

# Evoked controls (cluster split, glutamate dose, crosstalk): 60 runs, resumable
bash run/run_controls.sh 8
python3 -B scripts/plot_evoked_compare.py controls
```

See [model/README.md](model/README.md) for how the model files fit together and
[doc/README.md](doc/README.md) for details.
