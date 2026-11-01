"""Show free-calcium position along spine height throughout a 35 ms run.

The molecules include calcium from every modeled source and are not labeled
by NMDAR cluster. Color is the number in each 50 nm z bin, sampled every
0.1 ms and averaged over five frames to make the spatial pattern readable.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from plot_calcium_positions import positions
from cluster_cases import CASES, OUTPUT


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, blue_colormap, figure, save
apply_style()
Z_EDGES = np.arange(0.45, 1.751, 0.05)
ITERATIONS = np.arange(0, 350001, 1000)


def depth_counts(model_name: str, seed: int) -> np.ndarray:
    counts = []
    for iteration in ITERATIONS:
        calcium = positions(model_name, seed, int(iteration), 'Ca2')
        counts.append(np.histogram(calcium[:, 2], bins=Z_EDGES)[0])
    raw = np.asarray(counts, dtype=float)
    kernel = np.ones(5) / 5
    return np.apply_along_axis(lambda row: np.convolve(row, kernel, mode='same'),
                               axis=0, arr=raw)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--output-dir', type=Path,
                        default=OUTPUT)
    args = parser.parse_args()

    matrices = [depth_counts(model_name, args.seed) for model_name, _, _ in CASES]
    maximum = max(np.log10(1 + matrix).max() for matrix in matrices)
    for (model_name, _, _), matrix in zip(CASES, matrices):
        fig, ax = figure()
        image = ax.imshow(np.log10(1 + matrix.T), origin='lower', aspect='auto',
                          extent=(0, 35, Z_EDGES[0], Z_EDGES[-1]),
                          interpolation='nearest', cmap=blue_colormap(), vmin=0,
                          vmax=maximum)
        ax.set_ylabel('z (µm)')
        ax.set_xlabel('Time (ms)')
        bar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.04)
        bar.set_label('log10(1 + Ca ions / 50 nm)')
        save(fig, args.output_dir / 'conditions' / Path(model_name).name / 'ca_depth.png')


if __name__ == '__main__':
    main()
