# View the two calcium streams in Blender

The new source is `model/calcium_track/`: two clusters with 10 NMDARs
each, nominally 100 nm apart, without a spine apparatus. **Ca_A is pink** and
entered through cluster 1; **Ca_B is green** and entered through cluster 2.
Both have the same diffusion and chemistry as original calcium. **Ca2** is
background calcium, including VGCC entry and pump leaks.

## Corrected 200 nm follow-up

This run is **complete and validated**, including its clean two-species movie.

For the new 200 nm run, start a fresh scene and import:

```text
results/calcium_track/evoked/2_cluster_200nm/seed_00001/blender_ca_only/Scene.data_model.000000.json
```

In **Visualization Settings → Manually Select Viz Directory → Read Viz Data**,
select:

```text
results/calcium_track/evoked/2_cluster_200nm/seed_00001/blender_ca_only/frames/Scene.ascii.000000.dat
```

This movie has **3,501 frames**, spaced **0.01 ms (10 µs)** apart, over 35 ms.
Output index 1,000 is 10 ms; filename `Scene.ascii.100000.dat` is also 10 ms
because filenames encode MCell iterations. Do not use the 100 nm movie's
0.1 ms-per-frame conversion for this movie. The clean directory is generated
after simulation and validation finish; `status.json` reports completion.

The 200 nm cluster separation is mostly along **y**. Rotate the view or look
down onto the PSD to see both initial cluster locations; viewing only x–z
collapses most of their separation. Source colors remain pink/green, with only
free Ca_A/Ca_B particles in the clean movie.

Use the run ending **_002**. The earlier **_001** attempt was stopped because
the inherited fixed glutamate height put site 2 outside the local cleft. The
corrected release follows the curved membrane and passes the geometry check.
The existing 100 nm simulation remains unchanged and passes its placement audit.

## Open the clean Ca_A/Ca_B-only recording

A filtered copy of all 351 frames is ready. It contains only free **Ca_A** and
**Ca_B**, retaining the exact original particle IDs, coordinates and times.
No simulation was rerun. The full original recording remains available for
analysis. Future runs built from this model also generate this clean copy.

1. Open **Blender 4.4.3 CellBlender**, using a new scene for this model.
2. Choose **File → Import → CellBlender Model and Geometry (JSON)** and import:

   ```text
   results/calcium_track/evoked/2_cluster_100nm/seed_00001/blender_ca_only/Scene.data_model.000000.json
   ```

   This file includes geometry and the pink/green display settings. There should
   be no geometry object named `SA`. The `presynapse` mesh is a static boundary;
   the model does not simulate a presynaptic action potential.

3. In the 3D Viewport press **N** to open the sidebar. Open **CellBlender →
   Visualization Settings** and enable **Manually Select Viz Directory**.
4. Click **Read Viz Data** and select:

   ```text
   results/calcium_track/evoked/2_cluster_100nm/seed_00001/blender_ca_only/frames/Scene.ascii.000000.dat
   ```

   Select from the **frames/** subfolder, not its parent. This CellBlender
   version's manual loader enumerates every file in the chosen directory, so
   the JSON and manifest are deliberately kept outside the frames folder.
5. Play the timeline or select a frame. The complete test has **351 frames**
   including time zero. Output index 0–350 corresponds to **0–35 ms** in **0.1 ms**
   steps. For example, `Scene.ascii.100000.dat` is iteration 100,000 = **10 ms**.
   Filenames encode simulation iterations, not milliseconds.

Glutamate and the imposed voltage waveform begin together at t = 0. Labeled
calcium appears only after NMDARs conduct. The 1 ms startup test had no labeled
entry, so advance further into the full run to inspect both signals.

## Display settings and existing recordings

The MCell molecule **export/visualize checkboxes determine what is written during
a simulation**. Changing them after a run does not remove species from existing
particle files. For a reliable two-species movie, use the filtered recording
above rather than changing export settings on the original multi-species run.
Start from a new scene so previously loaded molecule objects do not remain.

- Show **Ca_A** and **Ca_B** in the molecule display controls. Imported RGB colors
  are pink `(1.00, 0.12, 0.55)` and green `(0.10, 0.90, 0.20)`. In a reused scene,
  check whether existing materials have overridden those colors.
- The clean recording contains **no Ca2, bound carriers or NMDAR particles**.
  No per-species hiding is needed. Their definitions remain in the import JSON
  to preserve the model, but their particle rows are absent from every frame.
- The original `viz_data/seed_00001/` recording includes background Ca2, bound
  carriers and all NMDAR states for analysis. Bound calcium is deliberately
  excluded from the clean movie because it is not freely propagating calcium.
- Use wireframe or a transparent material for the spine. Hide `Cube`, the static
  `presynapse`, and transparent counting surfaces if they obscure the calcium.
  Adjust molecule glyph size for visibility; glyphs are not physical ion radii.

Simultaneous glutamate release does not imply synchronized receptor opening or
continuous jets. The labels identify the **calcium entry cluster**, not the
release site of the glutamate that activated it. Calcium can bind, unbind,
diffuse down the neck into the dendrite, or be removed.

## Completion and spatial counts

From the repository root:

```bash
cat results/calcium_track/evoked/2_cluster_100nm/seed_00001/status.json
tail -n 20 results/calcium_track/evoked/2_cluster_100nm/seed_00001/simulation.log
tmux attach -t calcium-35ms
```

Detach with **Ctrl+b**, release, then **d**. The simulation runs under
`caffeinate -i`; keep the laptop powered with its lid open. The tmux session
exits when the job finishes. `state: complete` means simulation and automatic
validation/voxel analysis succeeded. For `state: failed`, inspect
`simulation.log` and `analysis.log`. Frames can be read during a running job,
but then represent only the completed portion.

The completed run contains:

```text
validation.json
spatial/calcium_voxels.csv.gz
spatial/metadata.json
```

The **50 nm voxels are postprocessing bins**, not barriers in the simulation.
They count free and bound calcium separately throughout the spine and attached
dendrite. No voxel overlay is automatically added to Blender. Original ASCII
files retain particle IDs and coordinates at every recorded instant. Events
between the 0.1 ms frames cannot be reconstructed exactly; see the model README
for 10 µs output settings and voxel coordinate conventions.

## Other runs and interpretation

To create a clean recording from another completed, validated run:

```bash
python3 -B model/calcium_track/export_clean_viz.py PATH_TO_RUN
```

This writes `PATH_TO_RUN/blender_ca_only/` and refuses to overwrite an existing
directory. `--seed N` selects another seed; `--destination PATH` chooses a new
output directory. The manifest verifies each frame's Ca_A/Ca_B counts against
the recorded count traces and records input/output hashes. Empty early frames
are retained so the timeline is not shifted.

For another run or seed, use both its own JSON and its own ASCII directory.
Do not combine a historical apparatus-containing JSON with these no-apparatus
particle files. Historical paths under `model/calcium_total/` remain intact;
their exports generally contain unlabelled `Ca2` and resting-state `NMDAR_0`.

“Evoked” and “spontaneous” remain provisional names for our glutamate+voltage
and glutamate-only protocols. Biological evoked release is triggered by a
presynaptic AP and does not require a postsynaptic backpropagating AP.
The [new model README](../model/README.md) records the
revised assumptions and deferred higher-glutamate tests.


## Overnight 200 nm repeats: seeds 2–5 (2026-10-15)

The batch is prepared at:
`results/calcium_track/evoked/2_cluster_200nm/`.
Check `batch_status.json` for completion. For each `seed_00002` through
`seed_00005`, use its `blender_ca_only/Scene.data_model.000000.json` and
`blender_ca_only/frames/` following the same import procedure above. These files
appear after that seed finishes validation and export. Each movie covers 35 ms
with 3,501 frames: frame 1000 is 10 ms. Pink is Ca_A; green is Ca_B.

`comparison/` will contain seed-1-through-5 calcium-entry and free-calcium plots,
a metrics CSV and a summary JSON. Each seed has its own `validation.json`,
`audit.json`, simulation log and source-resolved voxel counts.

Model scripts and input files now live outside results. `input_reference.json`
points to the matching frozen inputs. Existing seed-1 movies retain their paths.
The old input snapshots are archived; use the editable scripts under
`model/calcium_track/` for analysis and new runs.
