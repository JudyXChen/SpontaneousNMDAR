"""Voxel maps of where free Ca_A and Ca_B spend time, averaged across seeds.

Reads the voxel counts each calcium_track run already saved
(spatial/calcium_voxels.csv.gz); no simulation is run.

Each label is split by time since its own cluster last conducted. In each window:
  occupancy = free ions per voxel per frame, averaged across seeds (time windows use
              only seeds where that cluster opened)
  preference = (a - b) / (a + b), where a and b are each label's occupancy divided
               by its own total, so that a label entering more calcium does not
               look more widespread. -1 is all Ca_B (green), +1 is all Ca_A (pink),
               0 is an equal share. Voxels with few sampled ions are gray.

Maps are 2D: side view (the axis the clusters are separated along, against z, from the
spine root to the head) through a slab around the clusters, and top view (x, y) in a band just under the PSD. The ER
outline is drawn when the run has one. Voxels are clipped to the membrane outline only: the slab
also holds cytosol in front of and behind the ER, so color inside the ER line can be real.

Usage: plot_calcium_voxels.py OUTPUT_DIRECTORY RUN [RUN ...] [--state free|bound]
Bound maps include buffers, pumps and PLC, and read 'Bound ions' on the color bar.
"""
import argparse
import ast
import csv
import gzip
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path as PlotPath
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as style

ROOT = Path(__file__).resolve().parents[1]
PINK, GREEN = '#E52B88', '#169B46'
MASK = '#D9D9D9'
LABEL_MAP = LinearSegmentedColormap.from_list('ca_b_to_ca_a', (GREEN, '#FFFFFF', PINK))
MIN_SAMPLED_IONS = 20
# (name, lower ms, upper ms]: time since own cluster conducted; 'conducting' means exactly 0.
WINDOWS = (('conducting', 0.0, 0.0), ('0_to_50us', 0.0, 0.05),
           ('50us_to_1ms', 0.05, 1.0), ('after_1ms', 1.0, np.inf))
SLAB_HALF_WIDTH_UM = 0.15
BAND_BELOW_PSD_UM = (1.55, 1.70)
# preference = (r - 1) / (r + 1) for an enrichment ratio r, so the bar can read in fold change.
FOLD_TICKS = (-1, -3 / 5, -1 / 3, 0, 1 / 3, 3 / 5, 1)
FOLD_LABELS = ('only', '4×', '2×', 'equal', '2×', '4×', 'only')


def read_mesh(path, name):
    values = {node.targets[0].id: ast.literal_eval(node.value) for node in ast.parse(path.read_text()).body
              if isinstance(node, ast.Assign) and isinstance(node.value, ast.List)}
    return np.asarray(values[name + '_vertex_list']), np.asarray(values[name + '_wall_list'])


def mesh_section(vertices, walls, axis, value):
    """Line segments where a mesh crosses the plane coordinate[axis] == value."""
    segments = []
    for triangle in vertices[walls]:
        points = []
        for start, end in ((0, 1), (1, 2), (2, 0)):
            a, b = triangle[start], triangle[end]
            if (a[axis] - value) * (b[axis] - value) < 0:
                points.append(a + (value - a[axis]) / (b[axis] - a[axis]) * (b - a))
        if len(points) == 2:
            segments.append(points)
    return np.asarray(segments)


def outer_outline(segments):
    """Largest closed loop in a mesh section, as a clip path."""
    key = lambda point: tuple(np.round(point, 8))
    neighbors = {}
    for a, b in segments:
        neighbors.setdefault(key(a), []).append(key(b))
        neighbors.setdefault(key(b), []).append(key(a))
    loops, seen = [], set()
    for start in neighbors:
        if start in seen:
            continue
        loop, previous, current = [start], None, start
        seen.add(start)
        while True:
            following = [n for n in neighbors[current] if n != previous]
            if not following or following[0] == start or following[0] in seen:
                break
            previous, current = current, following[0]
            seen.add(current)
            loop.append(current)
        loops.append(np.asarray(loop))
    area = lambda xy: abs(np.dot(xy[:, 0], np.roll(xy[:, 1], 1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], 1)))
    best = max(loops, key=area)
    return PlotPath(np.vstack((best, best[:1])), closed=True)


def label_fold_bar(bar):
    bar.set_ticks(FOLD_TICKS, labels=FOLD_LABELS)
    bar.set_label('Enrichment, each label relative to its total')
    bar.ax.text(0, 1.2, 'Ca_B', transform=bar.ax.transAxes, color=GREEN, fontweight='bold', ha='left', va='bottom')
    bar.ax.text(1, 1.2, 'Ca_A', transform=bar.ax.transAxes, color=PINK, fontweight='bold', ha='right', va='bottom')


def time_since_conducting(run, seed, tag):
    trace = np.loadtxt(run / f'react_data/seed_{seed:05d}/NMDAR_{tag}_World.dat')
    times = trace[:, 0] * 1000
    since = np.full(len(times), np.inf)
    last = None
    for index, (time, conducting) in enumerate(zip(times, trace[:, 1])):
        if conducting > 0:
            last = time
        if last is not None:
            since[index] = time - last
    return since


def window_index(since):
    """Return the one time window a frame belongs to, as a tuple (empty before first opening)."""
    if not np.isfinite(since):
        return ()
    since = round(float(since), 6)
    if since == 0:
        return (0,)
    return tuple(i for i, (_, low, high) in enumerate(WINDOWS) if i and low < since <= high)


def load_seed(run, seed, along, slab_center, slab_half_width, state):
    """Voxel counts per window for one seed, as {(label, window, view): {(u, z): count}}."""
    since = {tag: time_since_conducting(run, seed, tag) for tag in 'AB'}
    frames = {tag: np.zeros(len(WINDOWS) + 1) for tag in 'AB'}
    for tag in 'AB':
        for value in since[tag]:
            frames[tag][-1] += 1
            for index in window_index(value):
                frames[tag][index] += 1
    counts = {}
    with gzip.open(run / 'spatial/calcium_voxels.csv.gz', 'rt', newline='') as stream:
        for row in csv.DictReader(stream):
            if row['state'] != state or row['source'] not in ('A', 'B'):
                continue
            tag = row['source']
            ix, iy, iz = int(row['ix']), int(row['iy']), int(row['iz'])
            count = int(row['count'])
            frame = int(row['iteration']) // 100
            windows = window_index(since[tag][frame]) + (len(WINDOWS),)
            size = 0.05
            index = (ix, iy)
            in_slab = abs((index[1 - along] + .5) * size - slab_center) <= slab_half_width
            in_band = BAND_BELOW_PSD_UM[0] <= (iz + .5) * size < BAND_BELOW_PSD_UM[1]
            for window in windows:
                if in_slab:
                    key = (tag, window, 'side', index[along], iz)
                    counts[key] = counts.get(key, 0) + count
                if in_band:
                    key = (tag, window, 'top', ix, iy)
                    counts[key] = counts.get(key, 0) + count
    return counts, frames


def accumulate(runs, along, slab_center, state):
    """Per-frame occupancy averaged across seeds, and raw sampled ions."""
    occupancy, sampled = {}, {}
    used = {(tag, w): 0 for tag in 'AB' for w in range(len(WINDOWS) + 1)}
    number = len(WINDOWS)
    for run in runs:
        seed = json.loads((run / 'status.json').read_text())['seed']
        print(f'Reading seed {seed}', flush=True)
        counts, frames = load_seed(run, seed, along, slab_center, SLAB_HALF_WIDTH_UM, state)
        for tag in 'AB':
            for window in range(number + 1):
                # A seed where a cluster never opened has no frames in its time windows.
                used[(tag, window)] += int(frames[tag][window] > 0)
        for (tag, window, view, i, j), count in counts.items():
            key = (tag, window, view, i, j)
            occupancy[key] = occupancy.get(key, 0.0) + count / frames[tag][window]
            sampled[key] = sampled.get(key, 0) + count
    return occupancy, sampled, used


def grids(values, tag, window, view, ranges):
    array = np.zeros((ranges[1][1] - ranges[1][0], ranges[0][1] - ranges[0][0]))
    for (t, w, v, i, j), value in values.items():
        if t == tag and w == window and v == view:
            if ranges[0][0] <= i < ranges[0][1] and ranges[1][0] <= j < ranges[1][1]:
                array[j - ranges[1][0], i - ranges[0][0]] = value
    return array


def save_map(fig, ax, image, path, colorbar_label=None):
    """Clean panel in fig/, viewing copy with colorbar, and a reusable colorbar PNG."""
    path = Path(path)
    (path.parent / 'fig').mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.parent / 'fig' / path.name, dpi=style.DPI)
    bar = fig.colorbar(image, ax=ax, orientation='horizontal', location='bottom', fraction=0.05, pad=0.2)
    if colorbar_label is None:
        label_fold_bar(bar)
    else:
        bar.set_label(colorbar_label)
    fig.savefig(path, dpi=style.DPI, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'Saved -> {path}')


def save_colorbar(image_map, norm, path):
    fig = plt.figure(figsize=(2.8, 0.6))
    ax = fig.add_axes((0.08, 0.62, 0.84, 0.22))
    label_fold_bar(fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=image_map), cax=ax, orientation='horizontal'))
    fig.savefig(path, dpi=style.DPI, transparent=True, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output', type=Path)
    parser.add_argument('runs', type=Path, nargs='+')
    parser.add_argument('--state', choices=('free', 'bound'), default='free')
    args = parser.parse_args()
    runs = [run.resolve() for run in args.runs]
    audit = json.loads((runs[0] / 'audit.json').read_text())
    centers = np.asarray(audit['region_centers_um'])
    A, B = centers
    separation = B[:2] - A[:2]
    along = int(abs(separation[1]) > abs(separation[0]))  # 0 = x, 1 = y
    axis_name = 'xy'[along]
    tilt = float(np.degrees(np.arctan2(abs(separation[1 - along]), abs(separation[along]))))
    slab_center = float((A[1 - along] + B[1 - along]) / 2)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    occupancy, sampled, used = accumulate(runs, along, slab_center, args.state)
    n = len(runs)
    size = 0.05
    vertices, walls = read_mesh(ROOT / 'model/shared/geometry/Mushroom_dendrite.py', 'Mushroom_dendrite')
    side_mesh = [mesh_section(vertices, walls, 1 - along, slab_center)]
    top_mesh = [mesh_section(vertices, walls, 2, 1.55)]
    clips = {'side': outer_outline(side_mesh[0][:, :, [along, 2]]), 'top': outer_outline(top_mesh[0][:, :, :2])}
    # Frozen inputs may stay on the computer that ran the simulation.
    if any('SA object present' in check for check in audit['checks']):
        er_vertices, er_walls = read_mesh(ROOT / 'model/shared/geometry/SA.py', 'SA')
        side_mesh.append(mesh_section(er_vertices, er_walls, 1 - along, slab_center))
        top_mesh.append(mesh_section(er_vertices, er_walls, 2, 1.5))
    views = {
        'side': dict(ranges=((-11, 10), (10, 35)),  # z from the spine root (0.5 um) to the head
                      xlabel=axis_name + ' (µm)', ylabel='z (µm)'),
        'top': dict(ranges=((-11, 10), (-11, 10)), xlabel='x (µm)', ylabel='y (µm)'),
    }
    style.apply_style()
    names = [w[0] for w in WINDOWS] + ['all_frames']
    summary = dict(seeds=n, state=args.state, voxel_um=size, side_axis=axis_name, tilt_of_separation_from_axis_deg=tilt, slab_center_um=slab_center,
                   slab_half_width_um=SLAB_HALF_WIDTH_UM, top_band_z_um=BAND_BELOW_PSD_UM,
                   min_sampled_ions=MIN_SAMPLED_IONS, windows={})
    arrays = {}
    for window, name in enumerate(names):
        summary['windows'][name] = {}
        for view, spec in views.items():
            ranges = spec['ranges']
            extent = (ranges[0][0] * size, ranges[0][1] * size, ranges[1][0] * size, ranges[1][1] * size)
            maps = {tag: grids(occupancy, tag, window, view, ranges) / max(used[(tag, window)], 1) for tag in 'AB'}
            raw = {tag: grids(sampled, tag, window, view, ranges) for tag in 'AB'}
            arrays[f'{view}_{name}_Ca_A'], arrays[f'{view}_{name}_Ca_B'] = maps['A'], maps['B']
            total = {tag: maps[tag].sum() for tag in 'AB'}
            enough = (raw['A'] + raw['B']) >= MIN_SAMPLED_IONS
            preference = np.full(maps['A'].shape, np.nan)
            if total['A'] > 0 and total['B'] > 0:
                a, b = maps['A'] / total['A'], maps['B'] / total['B']
                valid = enough & ((a + b) > 0)
                preference[valid] = (a[valid] - b[valid]) / (a[valid] + b[valid])
            arrays[f'{view}_{name}_preference'] = preference
            summary['windows'][name][view] = {
                'sampled_ions_A': int(raw['A'].sum()), 'sampled_ions_B': int(raw['B'].sum()),
                'voxels_shown': int(np.isfinite(preference).sum()),
                'seeds_averaged_A': used[('A', window)], 'seeds_averaged_B': used[('B', window)]}

            def draw(data, cmap, norm):
                fig, ax = style.figure()
                current = plt.get_cmap(cmap).copy() if isinstance(cmap, str) else cmap.copy()
                current.set_bad(MASK)
                image = ax.imshow(np.ma.masked_invalid(data), origin='lower', extent=extent, cmap=current,
                                  norm=norm, interpolation='nearest', aspect='equal')
                image.set_clip_path(PathPatch(clips[view], transform=ax.transData))
                for mesh in (side_mesh if view == 'side' else top_mesh):
                    for segment in mesh:
                        points = segment[:, [along, 2]] if view == 'side' else segment[:, :2]
                        ax.plot(*points.T, color='#333333', linewidth=0.5)
                ax.set(xlabel=spec['xlabel'], ylabel=spec['ylabel'], xlim=extent[:2], ylim=extent[2:])
                for center, color in ((A, PINK), (B, GREEN)):
                    point = (center[along], center[2]) if view == 'side' else (center[0], center[1])
                    ax.scatter(*point, marker='x', s=30, color='black', linewidths=1.8, zorder=4)
                    ax.scatter(*point, marker='x', s=22, color=color, linewidths=1.0, zorder=5)
                return fig, ax, image

            fig, ax, image = draw(preference, LABEL_MAP, Normalize(-1, 1))
            save_map(fig, ax, image, output / f'{view}_{name}_preference.png')
            top = max(maps['A'].max(), maps['B'].max())
            if top > 0:
                for tag in 'AB':
                    fig, ax, image = draw(maps[tag], style.blue_colormap(), Normalize(0, top))
                    save_map(fig, ax, image, output / f'{view}_{name}_Ca_{tag}.png',
                             f'{args.state.capitalize()} ions per voxel per frame')
    save_colorbar(LABEL_MAP, Normalize(-1, 1), output / 'fig/preference_colorbar.png')
    handles = [Line2D([], [], linestyle='', marker='x', color=PINK, label='Cluster A'),
               Line2D([], [], linestyle='', marker='x', color=GREEN, label='Cluster B'),
               Line2D([], [], color='#333333', linewidth=0.8, label='Membrane, ER')]
    style.save_legend(handles, output / 'legend.png', figsize=(3.2, 0.4), ncol=3)
    np.savez_compressed(output / 'voxel_maps.npz', **arrays)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary['windows'], indent=1))


if __name__ == '__main__':
    main()
