"""Plot matched one-cluster and two-cluster 35 ms MCell runs.

Counts are raw model outputs from one seed. Head/neck and base free-calcium
counts include every calcium source; only rxn_NMDAR firings identify modeled
NMDAR calcium entry. These curves do not separate the two cluster sources.
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
from scripts.plot_style import condition_color, apply_style, figure, save, save_legend
apply_style()
PANELS = (
    ('NMDAR_Ca_flux_World.dat', r'NMDAR Ca$^{2+}$ entry', 'ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_nmdar'),
    ('Ca2_Mushroom_transparent.dat',
     r'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_ca'),
    ('Ca2_head_neck.dat', 'Head and neck free calcium (ions)', 'ca_head_neck'),
    ('Ca2_Spine_base.dat', 'Spine-base free calcium (ions)', 'ca_base'),
)
CASE_COLORS = (
    condition_color('blue', 'bright'),
    condition_color('purple', 'bright'),
    condition_color('orange', 'bright'),
)


def plot_overlay(file_name: str, ylabel: str, suffix: str, seed: int,
                 series, output_dir: Path) -> None:
    fig, ax = figure()
    for (model_name, label, _), color in series:
        values = load_trace(model_name, seed, file_name)
        ax.plot(values[:, 0] * 1000, values[:, 1], color=color,
                linestyle='solid', label=label)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Time (ms)')
    ax.set_xlim(0, 35)
    save(fig, output_dir / f'{suffix}.png', legend_name=None)

    if suffix == 'ca_entry':
        scaled, scaled_ax = figure()
        for (model_name, label, _), color in series:
            values = load_trace(model_name, seed, file_name)
            scaled_ax.plot(values[:, 0] * 1000, values[:, 1],
                           color=color, linestyle='solid', label=label)
        scaled_ax.set_yscale('symlog', linthresh=10)
        scaled_ax.set_ylim(bottom=0)
        scaled_ax.set_ylabel(r'NMDAR Ca$^{2+}$ entry')
        scaled_ax.set_xlabel('Time (ms)')
        scaled_ax.set_xlim(0, 35)
        save(scaled, output_dir / 'ca_entry_symlog.png', legend_name=None)


def load_trace(model_name: str, seed: int, file_name: str) -> np.ndarray:
    path = results_for(model_for(model_name)) / 'react_data' / f'seed_{seed:05d}' / file_name
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2 or values[-1, 0] < 0.035 - 1e-8:
        raise ValueError(f'Expected a complete 35 ms two-column trace: {path}')
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--output-dir', type=Path,
                        default=OUTPUT)
    parser.add_argument('--all-comparison-dir', type=Path, default=None)
    parser.add_argument('--pair-comparison-dir', type=Path, default=None)
    args = parser.parse_args()
    all_comparison_dir = (
        args.all_comparison_dir.resolve() if args.all_comparison_dir
        else args.output_dir / ALL_COMPARISON
    )
    pair_comparison_dir = (
        args.pair_comparison_dir.resolve() if args.pair_comparison_dir
        else args.output_dir / PAIR_COMPARISON
    )

    all_series = tuple(zip(CASES, CASE_COLORS))
    # Both one-cluster controls use blue-family colors. The retained 10-NMDAR
    # one-cluster and 10+10-NMDAR two-cluster comparison uses the requested
    # original-kinetics blue and orange exactly.
    cluster_compare_series = (all_series[0], all_series[2])

    for file_name, ylabel, suffix in PANELS:
        for (model_name, label, _), color in all_series:
            values = load_trace(model_name, args.seed, file_name)
            single, single_ax = figure()
            single_ax.plot(values[:, 0] * 1000, values[:, 1], color=color,
                           linestyle='solid')
            single_ax.set_ylabel(ylabel)
            single_ax.set_xlabel('Time (ms)')
            single_ax.set_xlim(0, 35)
            save(single, args.output_dir / 'conditions' / Path(model_name).name
                 / f'{suffix}.png')
            print(f'{model_name} {file_name}: final={values[-1, 1]:g}, '
                  f'peak={values[:, 1].max():g}')

        plot_overlay(file_name, ylabel, suffix, args.seed, all_series,
                     all_comparison_dir)
        plot_overlay(file_name, ylabel, suffix, args.seed,
                     cluster_compare_series, pair_comparison_dir)

    handles = [Line2D([0], [0], color=color, linestyle='solid', label=label)
               for (_, label, _), color in all_series]
    save_legend(handles, all_comparison_dir / 'legend.png',
                figsize=(3.4, 0.9))
    compare_handles = [Line2D([0], [0], color=color, linestyle='solid', label=label)
                       for (_, label, _), color in cluster_compare_series]
    save_legend(compare_handles, pair_comparison_dir / 'legend.png',
                figsize=(3.4, 0.9))


if __name__ == '__main__':
    main()
