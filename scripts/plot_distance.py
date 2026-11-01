"""Compare calibrated 10+10 or 3+3 NMDAR runs at 50, 100, and 200 nm spacing."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

from matplotlib.lines import Line2D
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


PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'NMDAR_Ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_NMDAR'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'),
)


def load_trace(case: str, seed: int, file_name: str) -> np.ndarray:
    path = (ROOT / 'results/calcium_total/evoked/conditions' / case / 'react_data'
            / f'seed_{seed:05d}' / file_name)
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if len(values) != 3501 or not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(f'Expected a complete 35 ms trace with 3501 samples: {path}')
    return values


def save_panel(fig, stem: Path) -> None:
    save(fig, stem.with_suffix(".png"), legend_name=None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--receptors', type=int, choices=(3, 10), default=10)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    suffix = '_3_NMDAR_each' if args.receptors == 3 else ''
    cases = [(f'2_cluster_{d}nm{suffix}', f'{d} nm') for d in (50, 100, 200)]
    tone = 'bright' if args.receptors == 3 else 'dark'
    distance_colors = {case: condition_color(family, tone)
                       for (case, _), family in zip(cases, ('purple', 'orange', 'blue'))}
    comparison = 'compare_2cluster3NMDAR_dist' if args.receptors == 3 else 'compare_2cluster_dist'
    output_dir = (args.output_dir or ROOT / 'results/calcium_total/evoked' / comparison).resolve()
    apply_style()

    traces: dict[tuple[str, str], np.ndarray] = {}
    reference_time = None
    for file_name, _, suffix in PANELS:
        for case, _ in cases:
            values = load_trace(case, args.seed, file_name)
            if reference_time is None:
                reference_time = values[:, 0]
            elif not np.array_equal(reference_time, values[:, 0]):
                raise ValueError(f'Time grid differs for {case}/{file_name}')
            traces[(suffix, case)] = values

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / 'traces.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow([
            'time_s',
            *[f'{suffix}_{case}' for _, _, suffix in PANELS for case, _ in cases],
        ])
        for index, time_s in enumerate(reference_time):
            writer.writerow([
                f'{time_s:.8g}',
                *[f'{traces[(suffix, case)][index, 1]:.8g}'
                  for _, _, suffix in PANELS for case, _ in cases],
            ])
    print(f'Saved -> {csv_path}')

    for _, ylabel, suffix in PANELS:
        fig, ax = figure()
        for case, label in cases:
            values = traces[(suffix, case)]
            ax.plot(values[:, 0] * 1000, values[:, 1],
                    color=distance_colors[case],
                    drawstyle='steps-post', label=label)
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if suffix == 'open_NMDAR':
            ax.set_ylim(bottom=-0.05)
        save_panel(fig, output_dir / suffix)

    handles = [Line2D([], [], color=distance_colors[case], label=label)
               for case, label in cases]
    save_legend(handles, output_dir / 'legend.png', figsize=(3.2, 0.4), ncol=3)


if __name__ == '__main__':
    main()
