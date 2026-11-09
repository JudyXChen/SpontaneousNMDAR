"""Seed means with SD bands from every seed in results/calcium_total/{evoked,spont}/runs/CASE.

Figures:
  evoked/compare/1x10_vs_2x10/   1x10_500g vs 2x10_100nm_250g (EPSP + bAP)
  spont/compare/2x10_distance/   2x10 at 50/100/200 nm, 250 g (glutamate only, resting voltage)
The evoked distance series is plotted by plot_evoked_compare.py batch (evoked/compare/2x10_distance/).

Usage: plot_seed_means.py [--results DIR]
"""
import argparse
import csv
import json
from pathlib import Path
import sys

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, condition_color, figure, save, save_legend

TOTAL = ROOT / 'results/calcium_total'
PANELS = (('NMDAR_Ca_flux_World.dat', 'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'NMDAR_Ca_entry'),
          ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_NMDAR'),
          ('Ca2_Mushroom_transparent.dat', 'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'))
DISTANCE_COLORS = zip((50, 100, 200), ('purple', 'orange', 'blue'))


def find_seed_runs(results, protocol, case):
    """Return {seed: run directory} for every seed folder of one case."""
    return {int(work.name.split('_')[1]): work
            for work in sorted((results / protocol / 'runs' / case).glob('seed_[0-9][0-9][0-9][0-9][0-9]'))}


def load(run, seed, name, samples):
    path = next(run.rglob(f'react_data/seed_{seed:05d}/{name}'))
    values = np.loadtxt(path)
    if values.shape != (samples, 2) or not np.isfinite(values).all():
        raise ValueError(f'Incomplete trace {values.shape}, expected ({samples}, 2): {path}')
    return values


def make_figure(output, results, protocol, cases, samples, release_ms, xticks):
    """cases: (case folder, label, color). Writes one panel per observable, a legend and data."""
    output.mkdir(parents=True, exist_ok=True)
    found = {case: find_seed_runs(results, protocol, case) for case, _, _ in cases}
    if not all(len(runs) >= 1 for runs in found.values()):
        print(f'Skipping {output.name}; no finished seeds yet for ' +
              ', '.join(case for case, runs in found.items() if not runs))
        return
    summary = dict(protocol=protocol, std='sample SD across seeds, ddof=1; not SEM', cases={})
    rows = {}
    time_ms = None
    for name, ylabel, panel in PANELS:
        fig, ax = figure()
        for case, label, color in cases:
            traces = [load(run, seed, name, samples) for seed, run in found[case].items()]
            time_ms = traces[0][:, 0] * 1000 - release_ms
            stack = np.stack([trace[:, 1] for trace in traces])
            mean = stack.mean(axis=0)
            std = stack.std(axis=0, ddof=1) if len(traces) > 1 else np.zeros_like(mean)
            rows[f'{panel}_{case}_mean'], rows[f'{panel}_{case}_std'] = mean, std
            ax.fill_between(time_ms, np.maximum(0, mean - std), mean + std, color=color, alpha=0.16,
                            linewidth=0, step='post')
            ax.plot(time_ms, mean, color=color, drawstyle='steps-post')
            entry = summary['cases'].setdefault(case, dict(seeds=list(found[case]), final={}))
            entry['final'][panel] = dict(mean=float(mean[-1]), sample_sd=float(std[-1]),
                                         per_seed=dict(zip(map(str, found[case]), stack[:, -1].tolist())))
        xlabel = 'Time from glutamate release (ms)' if release_ms else 'Time (ms)'
        ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(time_ms[0], time_ms[-1]), ylim=(0, None))
        ax.set_xticks(xticks)
        save(fig, output / f'{panel}.png', legend_name=None)
    counts = sorted({len(runs) for runs in found.values()})
    n_label = ', '.join(map(str, counts))
    handles = [Line2D([], [], color=color, label=label) for _, label, color in cases]
    handles.append(Patch(facecolor='gray', alpha=0.16, edgecolor='none', label=f'±1 SD (n={n_label})'))
    save_legend(handles, output / 'legend.png', figsize=(3.2, 0.65), ncol=2)
    with (output / 'traces.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_ms', *rows])
        writer.writerows(zip(time_ms, *rows.values()))
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'{output.name}: ' + '; '.join(f'{case} seeds {list(runs)}' for case, runs in found.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--results', type=Path, default=TOTAL,
                        help='folder holding evoked/ and spont/; figures are written there too')
    args = parser.parse_args()
    results = args.results.resolve()
    apply_style()
    distances = [(f'2x10_{d}nm_250g', f'{d} nm', condition_color(family, 'dark')) for d, family in DISTANCE_COLORS]
    make_figure(results / 'evoked/compare/1x10_vs_2x10', results, 'evoked',
                [('1x10_500g', 'One cluster', condition_color('blue', 'dark')),
                 ('2x10_100nm_250g', 'Two clusters, 100 nm', condition_color('orange', 'dark'))],
                3501, 0, (0, 10, 20, 30))
    make_figure(results / 'spont/compare/2x10_distance', results, 'spont', distances, 4501, 10, (-10, 0, 10, 20, 30))


if __name__ == '__main__':
    main()
