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

MCell4 with its compatible Python, as shipped with [CellBlender](https://mcell.org).
`scripts/runtime.sh` finds a local Blender 4.4.3 CellBlender install automatically;
elsewhere (including WSL on Windows) set `MCELL_PATH` and `PYTHON`.

## Run

From the repository root:

```sh
# Evoked controls (cluster split, glutamate dose, crosstalk): 60 runs, resumable
bash run/run_controls.sh 8
python3 -B scripts/plot_evoked_compare.py controls   # also: batch, small-cluster
python3 -B scripts/plot_seed_means.py                # 1x10 vs 2x10; spont distance series

# One calcium_total seed by hand (evoked, then resting voltage)
python3 -B model/shared/run_legacy.py model/calcium_total/evoked/2x10_100nm_250g \
  results/calcium_total/evoked/runs/2x10_100nm_250g/seed_00006 6
python3 -B model/calcium_total/spont/2x10_100nm_250g/model.py --seeds 6 --output NEW_DIRECTORY

# Source-labeled calcium: build frozen inputs, then run one seed or a planned batch
python3 -B model/calcium_track/build.py model/calcium_track/snapshots/my_run \
  --condition model/calcium_track/conditions/2_cluster_200nm.json
python3 -B model/calcium_track/run.py results/calcium_track/my_run --inputs model/calcium_track/snapshots/my_run --seed 1
python3 -B model/calcium_track/batch.py model/calcium_track/batches/buffer_test_original_buffer.json --workers 5
```

Use the CellBlender Python (`$PYTHON`) for anything that imports MCell.
Earlier launchers (distance batches, buffer tests, tmux wrappers) are in Git history.

See [model/README.md](model/README.md) for how the model files fit together and
[doc/README.md](doc/README.md) for details.
