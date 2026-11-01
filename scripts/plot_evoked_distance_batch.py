"""Plot completed evoked distance batches as seed means with sample SD bands."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import sys

# Keep plotting and Python caches out of the repository.
os.environ.setdefault('MPLCONFIGDIR', '/tmp/spont-nmdar-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/spont-nmdar-cache')
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, condition_color, figure, save, save_legend

PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'NMDAR_Ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_NMDAR'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'),
)
DISTANCES = (50, 100, 200)
FAMILIES = ('purple', 'orange', 'blue')
GROUPS = (
    (10, '', 'dark', 'compare_2cluster_dist'),
    (3, '_3_NMDAR_each', 'bright', 'compare_2cluster3NMDAR_dist'),
)


def load_group(batch, seeds, receptor_count, case_suffix):
    expected_time = np.arange(3501) * 1e-5
    traces, warnings = {}, {}
    for distance in DISTANCES:
        case = f'2_cluster_{distance}nm{case_suffix}'
        for seed in seeds:
            work = batch / case / f'seed_{seed:05d}'
            for filename, _, panel in PANELS:
                path = work / 'react_data' / f'seed_{seed:05d}' / filename
                values = np.loadtxt(path)
                if values.shape != (len(expected_time), 2) or not np.isfinite(values).all():
                    raise ValueError(f'Invalid trace shape or nonfinite values: {path}')
                if not np.allclose(values[:, 0], expected_time, rtol=0, atol=1e-10):
                    raise ValueError(f'Incomplete or mismatched time grid: {path}')
                if np.any(values[:, 1] < 0):
                    raise ValueError(f'Negative recorded count: {path}')
                if panel == 'NMDAR_Ca_entry' and np.any(np.diff(values[:, 1]) < 0):
                    raise ValueError(f'Cumulative source counter decreases: {path}')
                if panel == 'open_NMDAR' and values[:, 1].max() > 2 * receptor_count:
                    raise ValueError(f'Open count exceeds receptor count: {path}')
                traces[(distance, seed, panel)] = values[:, 1]
            with (work / 'run.log').open('rb') as handle:
                handle.seek(max(0, (work / 'run.log').stat().st_size - 65536))
                ending = handle.read().decode(errors='replace')
            if 'Simulation CPU time =' not in ending:
                raise ValueError(f'Missing simulation completion report: {work / "run.log"}')
            warnings[f'{case}/seed_{seed:05d}'] = [
                line.strip() for line in ending.splitlines()
                if line.strip().startswith('Warning:') or '% of reactions missed.' in line
            ]
    return expected_time, traces, warnings


def plot_group(batch, output, seeds, receptor_count, tone, time, traces, warnings):
    output.mkdir(parents=True, exist_ok=True)
    colors = dict(zip(DISTANCES, [condition_color(family, tone) for family in FAMILIES]))
    stats = {}
    for _, _, panel in PANELS:
        for distance in DISTANCES:
            samples = np.stack([traces[(distance, seed, panel)] for seed in seeds])
            stats[(panel, distance)] = (samples.mean(axis=0), samples.std(axis=0, ddof=1))

    with (output / 'traces.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        columns = [(panel, distance, index) for _, _, panel in PANELS
                   for distance in DISTANCES for index in (0, 1)]
        writer.writerow(['time_s'] + [f'{panel}_{distance}nm_{("mean", "std")[index]}'
                                     for panel, distance, index in columns])
        for i, t in enumerate(time):
            writer.writerow([t] + [stats[(panel, distance)][index][i]
                                   for panel, distance, index in columns])

    handles = [Line2D([], [], color=colors[distance], label=f'{distance} nm')
               for distance in DISTANCES]
    handles.append(Patch(facecolor='gray', alpha=0.16, edgecolor='none',
                         label=f'±1 SD (n={len(seeds)})'))
    for _, ylabel, panel in PANELS:
        fig, ax = figure()
        for distance in DISTANCES:
            mean, std = stats[(panel, distance)]
            ax.fill_between(time * 1000, np.maximum(0, mean - std), mean + std,
                            color=colors[distance], alpha=0.16, linewidth=0,
                            step='post')
            ax.plot(time * 1000, mean, color=colors[distance],
                    drawstyle='steps-post', label=f'{distance} nm')
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, time[-1] * 1000), ylim=(0, None))
        ax.set_xticks((0, 10, 20, 30))
        ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.20),
                  ncol=2, frameon=False)
        save(fig, output / f'{panel}.png', legend_name=None)
    save_legend(handles, output / 'legend.png', figsize=(3.2, 0.65), ncol=2)

    final = {}
    for distance in DISTANCES:
        final[str(distance)] = {}
        for _, _, panel in PANELS:
            values = [float(traces[(distance, seed, panel)][-1]) for seed in seeds]
            final[str(distance)][panel] = {
                'per_seed_final': dict(zip(map(str, seeds), values)),
                'mean_final': float(np.mean(values)),
                'sample_std_final': float(np.std(values, ddof=1)),
            }
    metadata = {
        'batch': str(batch), 'seeds': seeds, 'n_seeds': len(seeds),
        'receptors_per_cluster': receptor_count, 'distances_nm': list(DISTANCES),
        'duration_s': float(time[-1]), 'samples_per_trace': len(time),
        'mean': 'arithmetic mean across the selected seeds at each recorded time',
        'std': 'sample standard deviation across seeds, ddof=1; not SEM',
        'band': 'mean ±1 sample SD, alpha=0.16; displayed lower bound clipped at zero',
        'time_series_stats': 'traces.csv retains unmodified means and sample SDs',
        'observables': {panel: filename for filename, _, panel in PANELS},
        'final_values': final, 'run_warning_summaries': warnings,
        'source_runs': {
            str(seed): 'historical seed 1, original placement order' if seed == 1
            else 'completed local batch, NMDAR-first placement order'
            for seed in seeds
        },
    }
    (output / 'summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for distance in DISTANCES:
        entry = final[str(distance)]['NMDAR_Ca_entry']
        print(f'{receptor_count}/cluster, {distance} nm: final entry '
              f'{entry["mean_final"]:.2f} ± {entry["sample_std_final"]:.2f} (sample SD)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=Path, help='result directory containing CASE/seed_NNNNN run folders')
    parser.add_argument('--seeds', type=int, nargs='+', default=[1, 2, 3, 4])
    parser.add_argument('--output-root', type=Path, default=ROOT / 'results/calcium_total/evoked')
    args = parser.parse_args()
    batch = args.batch.resolve()
    seeds = args.seeds
    if len(seeds) < 2 or len(set(seeds)) != len(seeds) or any(seed < 1 for seed in seeds):
        parser.error('Select at least two unique positive seeds')
    # Validate every input before exporting either comparison.
    loaded = [load_group(batch, seeds, count, suffix)
              for count, suffix, _, _ in GROUPS]
    apply_style()
    seed_label = '_seeds' + '_'.join(map(str, seeds))
    for (count, _, tone, folder), (time, traces, warnings) in zip(GROUPS, loaded):
        plot_group(batch, args.output_root.resolve() / (folder + seed_label),
                   seeds, count, tone, time, traces, warnings)


if __name__ == '__main__':
    main()
