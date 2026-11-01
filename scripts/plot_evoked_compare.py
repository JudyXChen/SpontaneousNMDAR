"""Compare calcium_total evoked conditions.

  plot_evoked_compare.py distance [--receptors 3|10] [--seed N]   seed-1 traces at 50/100/200 nm
  plot_evoked_compare.py batch BATCH_DIR [--seeds ...]            seed means with SD bands
  plot_evoked_compare.py small-cluster                            open NMDAR at low receptor counts
  plot_evoked_compare.py one-or-two                               one 10-NMDAR cluster vs 10+10
"""
from __future__ import annotations

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

EVOKED = ROOT / 'results/calcium_total/evoked'
CONDITIONS = EVOKED / 'conditions'
PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'NMDAR_Ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_NMDAR'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'),
)


def load_trace(path, require_samples=True):
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if (require_samples and len(values) != 3501) or not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(f'Expected a complete 35 ms trace: {path}')
    return values


def overlay_panels(cases, output_dir, seed, panels, legend_size, legend_columns, require_samples=True):
    """cases: (folder, label, color). Writes traces.csv, one panel per observable and a legend."""
    traces, reference = {}, None
    for file_name, _, panel in panels:
        for case, _, _ in cases:
            values = load_trace(CONDITIONS / case / 'react_data' / f'seed_{seed:05d}' / file_name, require_samples)
            if reference is None:
                reference = values[:, 0]
            elif not np.array_equal(reference, values[:, 0]):
                raise ValueError(f'Time grid differs for {case}/{file_name}')
            traces[(panel, case)] = values
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / 'traces.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_s', *[f'{panel}_{case}' for _, _, panel in panels for case, _, _ in cases]])
        for index, time_s in enumerate(reference):
            writer.writerow([f'{time_s:.8g}', *[f'{traces[(panel, case)][index, 1]:.8g}'
                                                for _, _, panel in panels for case, _, _ in cases]])
    for _, ylabel, panel in panels:
        fig, ax = figure()
        for case, label, color in cases:
            values = traces[(panel, case)]
            ax.plot(values[:, 0] * 1000, values[:, 1], color=color, drawstyle='steps-post', label=label)
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if panel == 'open_NMDAR':
            ax.set_ylim(bottom=-0.05)
        save(fig, output_dir / f'{panel}.png', legend_name=None)
    handles = [Line2D([], [], color=color, label=label) for _, label, color in cases]
    save_legend(handles, output_dir / 'legend.png', figsize=legend_size, ncol=legend_columns)


def run_distance(args):
    suffix = '_3_NMDAR_each' if args.receptors == 3 else ''
    tone = 'bright' if args.receptors == 3 else 'dark'
    cases = [(f'2_cluster_{d}nm{suffix}', f'{d} nm', condition_color(family, tone))
             for d, family in zip((50, 100, 200), ('purple', 'orange', 'blue'))]
    folder = 'compare_2cluster3NMDAR_dist' if args.receptors == 3 else 'compare_2cluster_dist'
    overlay_panels(cases, (args.output_dir or EVOKED / folder).resolve(), args.seed, PANELS, (3.2, 0.4), 3)


def run_one_or_two(args):
    panels = (('NMDAR_World.dat', 'Open NMDAR', 'open_NMDAR'),
              ('NMDAR_Ca_flux_World.dat', r'NMDAR Ca$^{2+}$ entry', 'NMDAR_Ca_entry'),
              ('Ca2_Mushroom_transparent.dat', r'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'))
    cases = [('1_cluster_10_NMDAR', 'One cluster - Calibrated kinetics', condition_color('blue', 'dark')),
             ('2_cluster_100nm', 'Two clusters - Calibrated kinetics', condition_color('orange', 'dark'))]
    overlay_panels(cases, (args.output_dir or EVOKED / 'compare_1or2cluster').resolve(), 1, panels, (4.2, 0.85), 2, False)


def run_small_cluster(args):
    output = (args.output_dir or EVOKED / 'compare_smallcluster').resolve()
    families = {1: 'blue', 3: 'purple', 5: 'orange'}
    for n_glu, tone in ((250, 'dark'), (500, 'bright')):
        fig, ax = figure()
        for n_nmdar, family in families.items():
            values = load_trace(CONDITIONS / f'1_cluster_{n_nmdar}_NMDAR_{n_glu}_glu/react_data/seed_00001/NMDAR_World.dat')
            ax.plot(values[:, 0] * 1000, values[:, 1], color=condition_color(family, tone),
                    drawstyle='steps-post', label=f'{n_nmdar} NMDAR')
        ax.set(xlabel='Time (ms)', ylabel='Open NMDAR (count)', xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        ax.set_ylim(bottom=-0.05)
        save(fig, output / f'open_NMDAR_{n_glu}_glutamate.png', legend_name=None)
        handles = [Line2D([], [], color=condition_color(family, tone), label=f'{count} NMDAR')
                   for count, family in families.items()]
        save_legend(handles, output / f'legend_{n_glu}_glutamate.png', figsize=(3.2, 0.4), ncol=3)


# ---- multi-seed batches ----

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


def run_batch(args):
    batch = args.batch.resolve()
    seeds = args.seeds
    if len(seeds) < 2 or len(set(seeds)) != len(seeds) or any(seed < 1 for seed in seeds):
        raise SystemExit('Select at least two unique positive seeds')
    # Validate every input before exporting either comparison.
    loaded = [load_group(batch, seeds, count, suffix) for count, suffix, _, _ in GROUPS]
    seed_label = '_seeds' + '_'.join(map(str, seeds))
    root = (args.output_dir or EVOKED).resolve()
    for (count, _, tone, folder), (time, traces, warnings) in zip(GROUPS, loaded):
        plot_group(batch, root / (folder + seed_label), seeds, count, tone, time, traces, warnings)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('distance')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--receptors', type=int, choices=(3, 10), default=10)
    p.set_defaults(function=run_distance)
    p.add_argument('--output-dir', type=Path)
    p = sub.add_parser('batch')
    p.add_argument('batch', type=Path, help='Directory containing CASE/seed_NNNNN run folders')
    p.add_argument('--seeds', type=int, nargs='+', default=[1, 2, 3, 4])
    p.add_argument('--output-dir', type=Path)
    p.set_defaults(function=run_batch)
    for name, function in (('small-cluster', run_small_cluster), ('one-or-two', run_one_or_two)):
        p = sub.add_parser(name)
        p.add_argument('--output-dir', type=Path)
        p.set_defaults(function=function)
    args = parser.parse_args()
    apply_style()
    args.function(args)


if __name__ == '__main__':
    main()
