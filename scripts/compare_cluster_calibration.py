"""Compare original and calibrated kinetics for one- and two-cluster runs.

Color encodes receptor arrangement:
  blue   = one cluster, 10 NMDAR
  orange = two clusters, 10 + 10 NMDAR

All traces use solid lines. Color shade distinguishes original from calibrated
kinetics within each cluster color family.

The script requires complete 35 ms seed-1 traces for every condition. It does
not modify or replace the earlier one-cluster or cluster-comparison figures.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.plot_style import (  # noqa: E402
    apply_style,
    condition_color,
    figure,
    save,
    save_legend,
)


SEED = 1
TRACE_DIRS = {
    ('one_cluster', 'original'):
        ROOT / 'archive/NMDAR_old/results/evoked/conditions/1_cluster_10_NMDAR/react_data/seed_00001',
    ('one_cluster', 'calibrated'):
        ROOT / 'results/calcium_total/evoked/conditions/1_cluster_10_NMDAR/react_data/seed_00001',
    ('two_cluster', 'original'):
        ROOT / 'archive/NMDAR_old/results/evoked/conditions/2_cluster_100nm/react_data/seed_00001',
    ('two_cluster', 'calibrated'):
        ROOT / 'results/calcium_total/evoked/conditions/2_cluster_100nm/react_data/seed_00001',
}
CLUSTERS = {
    'one_cluster': ('One cluster', 'blue'),
    'two_cluster': ('Two clusters', 'orange'),
}
KINETICS = {
    'original': ('Original kinetics', 'solid', 'main'),
    'calibrated': ('Calibrated kinetics', 'solid', 'dark'),
}
TRACE_COLORS = {
    ('one_cluster', 'original'): condition_color('blue', 'bright'),
    ('one_cluster', 'calibrated'): condition_color('blue', 'dark'),
    ('two_cluster', 'original'): condition_color('orange', 'bright'),
    ('two_cluster', 'calibrated'): condition_color('orange', 'dark'),
}
LINEWIDTHS = {'original': 1.6, 'calibrated': 1.6}
PANELS = {
    'open_NMDAR': ('NMDAR_World.dat', 'Open NMDAR'),
    'NMDAR_Ca_entry': ('NMDAR_Ca_flux_World.dat', r'NMDAR Ca$^{2+}$ entry'),
    # Free Ca2 inside the full spine cytosol. This includes all modeled
    # sources and sinks, unlike the cumulative NMDAR entry counter above.
    'spine_cytosolic_Ca': (
        'Ca2_Mushroom_transparent.dat',
        r'Spine cytosolic free Ca$^{2+}$ (ions)',
    ),
}


def load_trace(path: Path) -> np.ndarray:
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(
            f'Incomplete trace ends at {values[-1, 0]:g} s, expected 0.035 s: {path}')
    return values


def save_panel(fig, output_stem: Path) -> None:
    save(fig, output_stem.with_suffix(".png"), legend_name=None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        '--output-dir', type=Path, default=None,
        help='Override the hierarchy-aware results directory.')
    parser.add_argument(
        '--calibrated-only', action='store_true',
        help='Plot only the calibrated one- and two-cluster traces.')
    args = parser.parse_args()

    if args.output_dir is None:
        if args.calibrated_only:
            args.output_dir = (
                ROOT / 'results/calcium_total/evoked/'
                'compare_1or2cluster')
        else:
            args.output_dir = (
                ROOT / 'archive/NMDAR_old/comparison/'
                'evoked/compare_1_cluster_10_NMDAR_vs_2_cluster_100nm_seed1')

    selected_kinetics = (
        ('calibrated',) if args.calibrated_only
        else tuple(KINETICS)
    )

    apply_style()
    traces: dict[tuple[str, str, str], np.ndarray] = {}
    reference_time = None
    for panel, (file_name, _) in PANELS.items():
        for cluster in CLUSTERS:
            for kinetics in selected_kinetics:
                path = TRACE_DIRS[(cluster, kinetics)] / file_name
                values = load_trace(path)
                if reference_time is None:
                    reference_time = values[:, 0]
                elif not np.array_equal(values[:, 0], reference_time):
                    raise ValueError(f'Time grid differs from the other traces: {path}')
                traces[(panel, cluster, kinetics)] = values

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / 'traces.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow([
            'time_s',
            *[f'{panel}_{cluster}_{kinetics}'
              for panel in PANELS
              for cluster in CLUSTERS
              for kinetics in selected_kinetics],
        ])
        for index, time_s in enumerate(reference_time):
            writer.writerow([
                f'{time_s:.8g}',
                *[f'{traces[(panel, cluster, kinetics)][index, 1]:.8g}'
                  for panel in PANELS
                  for cluster in CLUSTERS
                  for kinetics in selected_kinetics],
            ])
    print(f'Saved -> {csv_path}')

    for panel, (_, ylabel) in PANELS.items():
        fig, ax = figure()
        for cluster in CLUSTERS:
            for kinetics in selected_kinetics:
                _, linestyle, _ = KINETICS[kinetics]
                values = traces[(panel, cluster, kinetics)]
                ax.plot(
                    values[:, 0] * 1e3,
                    values[:, 1],
                    color=TRACE_COLORS[(cluster, kinetics)],
                    linestyle=linestyle,
                    linewidth=LINEWIDTHS[kinetics],
                    drawstyle='steps-post',
                    label=f'{CLUSTERS[cluster][0]} - {KINETICS[kinetics][0]}',
                )
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if panel == 'open_NMDAR':
            ax.set_ylim(bottom=-0.05)
        save_panel(fig, args.output_dir / panel)

    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [], [],
            color=TRACE_COLORS[(cluster, kinetics)],
            linestyle=KINETICS[kinetics][1],
            linewidth=LINEWIDTHS[kinetics],
            label=f'{CLUSTERS[cluster][0]} - {KINETICS[kinetics][0]}',
        )
        for cluster in CLUSTERS
        for kinetics in selected_kinetics
    ]
    save_legend(
        handles,
        args.output_dir / 'legend.png',
        figsize=(4.2, 0.85),
        ncol=2,
    )


if __name__ == '__main__':
    main()
