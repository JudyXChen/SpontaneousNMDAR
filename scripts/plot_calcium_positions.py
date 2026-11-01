"""Show free-calcium positions along the spine at selected times.

The x-z projection does not identify which channel made each calcium ion.
NMDAR markers show receptor positions at time zero; later receptor states
may differ from NMDAR_0.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from cluster_cases import ALL_COMPARISON, CASES, OUTPUT, PAIR_COMPARISON


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for, model_for
from scripts.plot_style import apply_style, colors, figure, save, save_legend
apply_style()
TIMES_MS = (5, 15, 35)


def positions(model_name: str, seed: int, iteration: int, species: str) -> np.ndarray:
    path = (results_for(model_for(model_name)) / 'viz_data' / f'seed_{seed:05d}'
            / f'Scene.ascii.{iteration:06d}.dat')
    if iteration == 0 and not path.exists():
        path = path.with_name('Scene.ascii.0000.dat')
    points = []
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if fields[0] == species:
                points.append([float(value) for value in fields[2:5]])
    return np.asarray(points, dtype=float).reshape(-1, 3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--output-dir', type=Path,
                        default=OUTPUT)
    args = parser.parse_args()

    palette = colors(2)
    for model_name, _, _ in CASES:
        receptors = positions(model_name, args.seed, 0, 'NMDAR_0')
        for time_ms in TIMES_MS:
            calcium = positions(model_name, args.seed, time_ms * 10000, 'Ca2')
            fig, ax = figure()
            ax.scatter(calcium[:, 0], calcium[:, 2], s=3, alpha=0.45,
                       color=palette[0], label='Free Ca')
            ax.scatter(receptors[:, 0], receptors[:, 2], s=12, alpha=0.8,
                       color=palette[1], label='NMDAR at 0 ms')
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(0.45, 1.75)
            ax.set_xlabel('x (µm)')
            ax.set_ylabel('z (µm)')
            print(f'{model_name}, {time_ms} ms: {len(calcium)} free Ca ions')
            save(fig, args.output_dir / 'conditions' / Path(model_name).name
                 / f'ca_{time_ms}ms.png', legend_name='positions_legend.png')

    handles = [Line2D([0], [0], linestyle='none', marker='o', markersize=5,
                      markerfacecolor=color, markeredgecolor=color,
                      label=label)
               for color, label in zip(palette, ('Free Ca', 'NMDAR at 0 ms'))]
    save_legend(handles, args.output_dir / ALL_COMPARISON / 'positions_legend.png',
                figsize=(3, 0.35), ncol=2, handlelength=1)
    save_legend(handles, args.output_dir / PAIR_COMPARISON / 'positions_legend.png',
                figsize=(3, 0.35), ncol=2, handlelength=1)


if __name__ == '__main__':
    main()
