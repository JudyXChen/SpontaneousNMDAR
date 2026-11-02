"""Side-view scatter of free Ca_A and Ca_B at one time, one figure per seed.

Same panel as compare_spatial.py's preview, for every seed_NNNNN folder in a batch.

Runs that kept ASCII frames are drawn from their exact coordinates. Runs that kept only
spatial/calcium_voxels.csv.gz know each ion's position to 50 nm, so each ion is drawn at a
random point in its voxel, restricted to cytosol (inside the spine and dendrite mesh,
outside the ER). Voxel counts are exact; positions within a voxel are not data. The
`positions` column of counts.csv says which applies, and the random draw is fixed by seed.

Usage: plot_calcium_positions.py BATCH [--time-ms 10] [--output BATCH/calcium_positions]
"""
import argparse
import ast
import csv
import gzip
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as style

GEOMETRY = Path(__file__).resolve().parents[1] / 'model/shared/geometry'
COLORS = {'A': '#E52B88', 'B': '#169B46'}
TIME_STEP_S = 1e-7
VOXEL_UM = 0.05


def read_mesh(name):
    values = {node.targets[0].id: ast.literal_eval(node.value)
              for node in ast.parse((GEOMETRY / f'{name}.py').read_text()).body
              if isinstance(node, ast.Assign) and isinstance(node.value, ast.List)}
    return np.asarray(values[name + '_vertex_list']), np.asarray(values[name + '_wall_list'])


def inside(points, mesh):
    """Parity test along +z. Triangles are found by xy barycentric coordinates."""
    triangles = mesh[0][mesh[1]]
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    keep = np.abs(det) > 1e-14
    a, b, c, det = a[keep], b[keep], c[keep], det[keep]
    px, py = points[:, 0, None], points[:, 1, None]
    u = ((px - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (py - a[:, 1])) / det
    v = ((b[:, 0] - a[:, 0]) * (py - a[:, 1]) - (px - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
    z = a[:, 2] + u * (b[:, 2] - a[:, 2]) + v * (c[:, 2] - a[:, 2])
    hits = (u >= 0) & (v >= 0) & (u + v <= 1) & (z > points[:, 2, None])
    return hits.sum(axis=1) % 2 == 1


def exact_positions(run, seed, iteration):
    """Coordinates by label from a saved ASCII frame, or None if the run kept none."""
    for folder in (run / f'viz_data/seed_{seed:05d}', run / 'blender_ca_only/frames'):
        path = folder / f'Scene.ascii.{iteration:06d}.dat'
        if path.exists():
            rows = {'A': [], 'B': []}
            for line in path.read_text().splitlines():
                fields = line.split()
                if fields and fields[0] in ('Ca_A', 'Ca_B'):
                    rows[fields[0][-1]].append([float(value) for value in fields[2:5]])
            return {tag: np.asarray(values).reshape(-1, 3) for tag, values in rows.items()}
    return None


def voxel_counts(run, iteration):
    """{(label, ix, iy, iz): free ions} at one saved frame."""
    counts = {}
    with gzip.open(run / 'spatial/calcium_voxels.csv.gz', 'rt', newline='') as stream:
        for row in csv.DictReader(stream):
            frame = int(row['iteration'])
            if frame > iteration:
                break
            if frame == iteration and row['state'] == 'free' and row['source'] in ('A', 'B'):
                counts[(row['source'], int(row['ix']), int(row['iy']), int(row['iz']))] = int(row['count'])
    return counts


def resample(counts, rng, cytosol, er):
    """Draw each ion uniformly over the cytosol in its voxel; returns positions, unplaced count."""
    rows = {'A': [], 'B': []}
    unplaced = 0
    for (tag, *index), count in counts.items():
        corner = np.asarray(index) * VOXEL_UM
        candidates = corner + rng.random((max(400, 40 * count), 3)) * VOXEL_UM
        candidates = candidates[inside(candidates, cytosol) & ~inside(candidates, er)]
        if len(candidates) < count:  # voxel barely touches cytosol
            unplaced += count - len(candidates)
            candidates = np.vstack((candidates, corner + rng.random((count - len(candidates), 3)) * VOXEL_UM))
        rows[tag].append(candidates[:count])
    return {tag: np.vstack(values) if values else np.empty((0, 3)) for tag, values in rows.items()}, unplaced


def draw(xyz, centers, wireframe, path):
    mid = centers.mean(axis=0)
    horizontal = (centers[1] - centers[0])[:2]
    horizontal /= np.linalg.norm(horizontal)
    vertices, walls = wireframe
    projected = np.column_stack(((vertices[:, :2] - mid[:2]) @ horizontal, vertices[:, 2]))
    edges = np.concatenate([walls[:, [0, 1]], walls[:, [1, 2]], walls[:, [2, 0]]])
    fig, ax = style.figure()
    ax.add_collection(LineCollection(projected[edges], color='#999999', linewidths=.2, alpha=.15))
    for tag, points in xyz.items():
        ax.scatter((points[:, :2] - mid[:2]) @ horizontal, points[:, 2], s=7, color=COLORS[tag], alpha=.65,
                   linewidths=0, label='Ca_' + tag)
    for tag, center in zip('AB', centers):
        ax.scatter([(center[:2] - mid[:2]) @ horizontal], [center[2]], marker='x', s=35, color=COLORS[tag],
                   linewidths=1.3)
    ax.set(xlabel='Along cluster separation (µm)', ylabel='z (µm)', xlim=(-.75, .75), ylim=(.35, 1.75))
    ax.set_aspect('equal')
    style.save(fig, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('batch', type=Path)
    parser.add_argument('--time-ms', type=float, default=10.0)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = args.output or args.batch / 'calcium_positions'
    iteration = round(args.time_ms * 1e-3 / TIME_STEP_S)
    runs = sorted(path for path in args.batch.glob('seed_*') if path.is_dir())
    if not runs:
        raise ValueError(f'No seed folders in {args.batch}')
    style.apply_style()
    wireframe, cytosol, er = read_mesh('Mushroom_transparent'), read_mesh('Mushroom_dendrite'), read_mesh('SA')
    table = []
    for run in runs:
        seed = json.loads((run / 'status.json').read_text())['seed']
        centers = np.asarray(json.loads((run / 'audit.json').read_text())['region_centers_um'])
        xyz, unplaced, source = exact_positions(run, seed, iteration), 0, 'exact'
        if xyz is None:
            xyz, unplaced = resample(voxel_counts(run, iteration), np.random.default_rng(seed), cytosol, er)
            source = 'random within 50 nm voxel'
        draw(xyz, centers, wireframe, output / f'seed_{seed:05d}_{args.time_ms:.3f}ms.png')
        table.append(dict(seed=seed, time_ms=args.time_ms, Ca_A=len(xyz['A']), Ca_B=len(xyz['B']),
                          positions=source, placed_outside_cytosol=unplaced))
        print(table[-1], flush=True)
    with (output / 'counts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)


if __name__ == '__main__':
    main()
