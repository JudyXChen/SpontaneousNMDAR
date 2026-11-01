"""Compare open NMDAR counts for one-cluster low-count dose runs."""

from __future__ import annotations

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


OUTPUT = ROOT / 'results/calcium_total/evoked/compare_smallcluster'
NMDAR_COUNTS = (1, 3, 5)
GLUTAMATE_COUNTS = (250, 500)
COUNT_FAMILIES = {1: 'blue', 3: 'purple', 5: 'orange'}


def load_trace(n_nmdar: int, n_glu: int) -> np.ndarray:
    path = (
        ROOT / 'results/calcium_total/evoked/conditions'
        / f'1_cluster_{n_nmdar}_NMDAR_{n_glu}_glu'
        / 'react_data/seed_00001/NMDAR_World.dat'
    )
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if len(values) != 3501 or not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(f'Expected a complete 35 ms trace with 3501 samples: {path}')
    return values


def main() -> None:
    apply_style()
    for n_glu in GLUTAMATE_COUNTS:
        tone = 'dark' if n_glu == 250 else 'bright'
        fig, ax = figure()
        for n_nmdar in NMDAR_COUNTS:
            values = load_trace(n_nmdar, n_glu)
            color = condition_color(COUNT_FAMILIES[n_nmdar], tone)
            ax.plot(values[:, 0] * 1000, values[:, 1], color=color,
                    drawstyle='steps-post', label=f'{n_nmdar} NMDAR')
        ax.set(xlabel='Time (ms)', ylabel='Open NMDAR (count)', xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        ax.set_ylim(bottom=-0.05)
        save(fig, OUTPUT / f'open_NMDAR_{n_glu}_glutamate.png', legend_name=None)

    for n_glu in GLUTAMATE_COUNTS:
        tone = 'dark' if n_glu == 250 else 'bright'
        handles = [
            Line2D([], [], color=condition_color(COUNT_FAMILIES[count], tone),
                   label=f'{count} NMDAR')
            for count in NMDAR_COUNTS
        ]
        save_legend(handles, OUTPUT / f'legend_{n_glu}_glutamate.png',
                    figsize=(3.2, 0.4), ncol=3)


if __name__ == '__main__':
    main()
