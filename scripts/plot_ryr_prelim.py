"""Per-seed RyR opening and ER calcium in the clamped buffer test. No simulation is run.

Seeds are drawn individually because RyR opening happens in only some of them, so a mean
would describe no actual run.

Usage: plot_ryr_prelim.py [--output DIR]
"""
import argparse
import json
from pathlib import Path
import sys

from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, figure, save, save_legend, shade

TEST = ROOT / 'results/calcium_track/clamped/buffer_test_200nm_1000glu'
BUFFERS = (('original_buffer', 'Membrane buffer (original)', shade('blue', 'main')),
           ('cytosolic_buffer', 'Cytosolic buffer', shade('purple', 'main')))
SEEDS = range(1, 6)


def trace(buffer, seed, name):
    return np.loadtxt(TEST / buffer / f'seed_{seed:05d}/react_data/seed_{seed:05d}/{name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=TEST / 'ryr_prelim')
    args = parser.parse_args()
    apply_style()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    panels = (('open_RyR', 'Open RyR (count)', lambda b, s: trace(b, s, 'RyR_O_World.dat')),
              ('ER_Ca', 'ER Ca$^{2+}$ (ions)', lambda b, s: trace(b, s, 'Ca2_ER_SA.dat')),
              ('labeled_Ca_to_RyR', 'Cumulative labeled Ca$^{2+}$\nbound to RyR (ions)',
               lambda b, s: np.column_stack((trace(b, s, 'Ca_A_to_RyR_World.dat')[:, 0],
                                             trace(b, s, 'Ca_A_to_RyR_World.dat')[:, 1]
                                             + trace(b, s, 'Ca_B_to_RyR_World.dat')[:, 1]))))
    for panel, ylabel, read in panels:
        fig, ax = figure()
        for buffer, _, color in BUFFERS:
            for seed in SEEDS:
                values = read(buffer, seed)
                ax.plot(values[:, 0] * 1000, values[:, 1], color=color, drawstyle='steps-post', linewidth=1.2)
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35), ylim=(0, None))
        ax.set_xticks((0, 10, 20, 30))
        save(fig, output / f'{panel}.png', legend_name=None)
    handles = [Line2D([], [], color=color, label=f'{label} (5 seeds)') for _, label, color in BUFFERS]
    save_legend(handles, output / 'legend.png', figsize=(3.2, 0.65), ncol=1)
    summary = {}
    for buffer, _, _ in BUFFERS:
        summary[buffer] = {}
        for seed in SEEDS:
            opened = trace(buffer, seed, 'RyR_O_World.dat')
            er = trace(buffer, seed, 'Ca2_ER_SA.dat')
            entry = sum(trace(buffer, seed, f'NMDAR_{tag}_Ca_flux_World.dat')[-1, 1] for tag in 'AB')
            summary[buffer][seed] = dict(
                first_open_ms=float(opened[np.argmax(opened[:, 1] > 0), 0] * 1000) if opened[:, 1].max() > 0 else None,
                open_at_35ms=float(opened[-1, 1]), entry_at_35ms=float(entry),
                er_start=float(er[0, 1]), er_at_35ms=float(er[-1, 1]))
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
