"""Plot the realized seed-1 NMDAR positions and glutamate sites on the PSD."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from plot_calcium_positions import positions
from cluster_cases import CASES, OUTPUT


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, colors, figure, legend_below, save
apply_style()


def mesh_data() -> tuple[np.ndarray, list[list[int]], list[int]]:
    path = (ROOT / 'archive/NMDAR_old/model/evoked/2_cluster_100nm/geometry/'
            'Mushroom_dendrite.py')
    names = {'Mushroom_dendrite_vertex_list', 'Mushroom_dendrite_wall_list',
             'Mushroom_dendrite_PSD_wall_indices'}
    values = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in names:
                values[target.id] = ast.literal_eval(node.value)
    if values.keys() != names:
        raise ValueError(f'Missing PSD geometry in {path}: {names - values.keys()}')
    return (np.asarray(values['Mushroom_dendrite_vertex_list']),
            values['Mushroom_dendrite_wall_list'],
            values['Mushroom_dendrite_PSD_wall_indices'])


def polygons(vertices: np.ndarray, walls: list[list[int]], indices) -> np.ndarray:
    return np.asarray([vertices[walls[index], :2] * 1000 for index in indices])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--output-dir', type=Path,
                        default=OUTPUT)
    args = parser.parse_args()

    vertices, walls, psd_indices = mesh_data()
    base = polygons(vertices, walls, psd_indices)
    palette = colors(3)
    for model_name, label, regions in CASES:
        if isinstance(regions[0], int):
            regions = (regions,)
        centers = [polygons(vertices, walls, indices).mean(axis=(0, 1))
                   for indices in regions]
        fig, ax = figure()
        ax.add_collection(PolyCollection(base, facecolor='#eeeeee',
                                         edgecolor='#d2d2d2', linewidth=0.3))
        receptors = positions(model_name, args.seed, 0, 'NMDAR_0')[:, :2] * 1000
        for number, (indices, center, color) in enumerate(
                zip(regions, centers, palette[:2]), 1):
            ax.add_collection(PolyCollection(polygons(vertices, walls, indices),
                                             facecolor=color, alpha=0.2,
                                             edgecolor=color, linewidth=1))
            assigned = np.argmin(np.linalg.norm(
                receptors[:, None, :] - np.asarray(centers)[None, :, :], axis=2),
                axis=1) == number - 1
            ax.scatter(receptors[assigned, 0], receptors[assigned, 1],
                       s=24, color=color, label=f'Cluster {number} (n={assigned.sum()})',
                       zorder=3)
            ax.scatter(center[0], center[1], s=90, marker='x', c=[palette[2]],
                       linewidths=2, label='Glutamate site' if number == 1 else None,
                       zorder=4)
        if len(centers) == 2:
            print(f'Cluster-center separation: {np.linalg.norm(centers[0] - centers[1]):.1f} nm')
        ax.set_aspect('equal')
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')
        legend_below(ax)
        save(fig, args.output_dir / 'conditions' / Path(model_name).name / 'psd_layout.png')


if __name__ == '__main__':
    main()
