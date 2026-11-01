"""Plot one completed 35 ms evoked NMDAR-cluster case without comparisons."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for, ARCHIVE
from scripts.plot_style import (  # noqa: E402
    condition_output_dir,
    apply_style,
    colors,
    condition_color,
    figure,
    legend_below,
    save,
)


PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_nmdar'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_ca'),
    ('Ca2_head_neck.dat', 'Head and neck free Ca$^{2+}$ (ions)', 'ca_head_neck'),
    ('Ca2_Spine_base.dat', 'Spine-base free Ca$^{2+}$ (ions)', 'ca_base'),
)

TRACE_COLORS = {
    'NMDAR_old/evoked/1_cluster_10_NMDAR': condition_color('blue', 'bright'),
    'NMDAR_old/evoked/1_cluster_20_NMDAR': condition_color('purple', 'bright'),
    'NMDAR_old/evoked/2_cluster_100nm': condition_color('orange', 'bright'),
    'calcium_total/evoked/1_cluster_10_NMDAR': condition_color('blue', 'dark'),
    'calcium_total/evoked/2_cluster_50nm': condition_color('purple', 'dark'),
    'calcium_total/evoked/2_cluster_100nm': condition_color('orange', 'dark'),
    'calcium_total/evoked/2_cluster_200nm': condition_color('blue', 'dark'),
    'calcium_total/evoked/2_cluster_50nm_3_NMDAR_each':
        condition_color('purple', 'dark'),
    'calcium_total/evoked/2_cluster_100nm_3_NMDAR_each':
        condition_color('orange', 'dark'),
    'calcium_total/evoked/2_cluster_200nm_3_NMDAR_each':
        condition_color('blue', 'dark'),
}


def trace_color(model_dir: Path) -> str:
    """Return the approved condition color, or the default for other cases."""
    try:
        key = model_dir.relative_to(ROOT / 'model').as_posix()
    except ValueError:
        try:
            key = 'NMDAR_old/' + model_dir.relative_to(ARCHIVE / 'model').as_posix()
        except ValueError:
            return colors(1)[0]
    return TRACE_COLORS.get(key, colors(1)[0])


def literal_assignments(path: Path, names: set[str]) -> dict[str, object]:
    values: dict[str, object] = {}
    for node in ast.parse(path.read_text()).body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in names):
            values[node.targets[0].id] = ast.literal_eval(node.value)
    missing = names - values.keys()
    if missing:
        raise ValueError(f'Missing assignments in {path}: {sorted(missing)}')
    return values


def load_trace(model_dir: Path, seed: int, file_name: str) -> np.ndarray:
    path = results_for(model_dir) / 'react_data' / f'seed_{seed:05d}' / file_name
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if len(values) != 3501 or not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(f'Expected a complete 35 ms trace with 3501 samples: {path}')
    return values


def positions(model_dir: Path, seed: int, iteration: int, species: str) -> np.ndarray:
    path = (results_for(model_dir) / 'viz_data' / f'seed_{seed:05d}'
            / f'Scene.ascii.{iteration:06d}.dat')
    if iteration == 0 and not path.exists():
        path = path.with_name('Scene.ascii.0000.dat')
    points = []
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0] == species:
                points.append([float(value) for value in fields[2:5]])
    return np.asarray(points, dtype=float).reshape(-1, 3)


def plot_layout(model_dir: Path, output_dir: Path, seed: int) -> None:
    geometry = literal_assignments(
        model_dir / 'geometry/Mushroom_dendrite.py',
        {'Mushroom_dendrite_vertex_list', 'Mushroom_dendrite_wall_list',
         'Mushroom_dendrite_PSD_wall_indices'},
    )
    layout = literal_assignments(
        model_dir / 'instantiation.py', {'CLUSTER_WALL_INDICES'})
    vertices = np.asarray(geometry['Mushroom_dendrite_vertex_list'])
    walls = geometry['Mushroom_dendrite_wall_list']
    regions = layout['CLUSTER_WALL_INDICES']
    if regions and isinstance(regions[0], int):
        regions = (regions,)

    def polygons(indices) -> np.ndarray:
        return np.asarray([vertices[walls[index], :2] * 1000 for index in indices])

    centers = np.asarray([polygons(indices).mean(axis=(0, 1)) for indices in regions])
    receptors = positions(model_dir, seed, 0, 'NMDAR_0')[:, :2] * 1000
    assigned = np.argmin(
        np.linalg.norm(receptors[:, None, :] - centers[None, :, :], axis=2), axis=1)
    palette = colors(3)
    fig, ax = figure()
    ax.add_collection(PolyCollection(
        polygons(geometry['Mushroom_dendrite_PSD_wall_indices']),
        facecolor='#eeeeee', edgecolor='#d2d2d2', linewidth=0.3))
    for number, (indices, center, color) in enumerate(
            zip(regions, centers, palette[:2]), 1):
        selected = assigned == number - 1
        ax.add_collection(PolyCollection(
            polygons(indices), facecolor=color, alpha=0.2,
            edgecolor=color, linewidth=1))
        ax.scatter(receptors[selected, 0], receptors[selected, 1], s=24,
                   color=color, label=f'Cluster {number} (n={selected.sum()})',
                   zorder=3)
        ax.scatter(center[0], center[1], s=90, marker='x', color=palette[2],
                   linewidths=2,
                   label='Glutamate site' if number == 1 else None, zorder=4)
    if len(centers) == 2:
        separation = np.linalg.norm(centers[0] - centers[1])
        print(f'Realized cluster-center separation: {separation:.3f} nm')
    ax.set_aspect('equal')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    legend_below(ax)
    save(fig, output_dir / 'psd_layout.png')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model_dir', type=Path)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--output-dir', type=Path, default=None)
    args = parser.parse_args()
    model_dir = args.model_dir.resolve()
    if args.output_dir is None:
        output_dir = condition_output_dir(model_dir)
    else:
        output_dir = args.output_dir.resolve()

    apply_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_layout(model_dir, output_dir, args.seed)
    color = trace_color(model_dir)
    for file_name, ylabel, suffix in PANELS:
        values = load_trace(model_dir, args.seed, file_name)
        fig, ax = figure()
        ax.plot(values[:, 0] * 1000, values[:, 1], color=color,
                drawstyle='steps-post')
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if file_name == 'NMDAR_World.dat':
            ax.set_ylim(bottom=-0.05)
        save(fig, output_dir / f'{suffix}.png')
        print(f'{file_name}: final={values[-1, 1]:g}, peak={values[:, 1].max():g}')


if __name__ == '__main__':
    main()
