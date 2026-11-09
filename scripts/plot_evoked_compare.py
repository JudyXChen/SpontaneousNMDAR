"""Compare calcium_total evoked conditions.

  plot_evoked_compare.py batch [--seeds ...]                 2x10 and 2x3 distance series, mean ± SD over seeds
  plot_evoked_compare.py controls                            run_controls.sh conditions, seeds 1-5, with open RyR
  plot_evoked_compare.py distance [--receptors 3|10] [--seed N]   one seed at 50/100/200 nm
  plot_evoked_compare.py one-or-two [--seed N]               1x10_500g vs 2x10_100nm_250g, one seed
  plot_evoked_compare.py small-cluster                       open NMDAR for 1x1, 1x3, 1x5, seed 1

Case names: <clusters>x<NMDAR per cluster>[_<distance>nm]_<glutamate per site>g[_Aonly].
Reads results/calcium_total/evoked/runs/CASE/seed_NNNNN; writes results/calcium_total/evoked/compare/NAME.
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
from scripts.plot_style import apply_style, colors, condition_color, figure, save, save_legend

EVOKED = ROOT / 'results/calcium_total/evoked'
RUNS = EVOKED / 'runs'
COMPARE = EVOKED / 'compare'
PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'NMDAR_Ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_NMDAR'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'),
)


def trace_path(case, seed, file_name, runs=RUNS):
    return runs / case / f'seed_{seed:05d}' / 'react_data' / f'seed_{seed:05d}' / file_name


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
            values = load_trace(trace_path(case, seed, file_name), require_samples)
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
    tone = 'bright' if args.receptors == 3 else 'dark'
    cases = [(f'2x{args.receptors}_{d}nm_250g', f'{d} nm', condition_color(family, tone))
             for d, family in zip((50, 100, 200), ('purple', 'orange', 'blue'))]
    folder = f'2x{args.receptors}_distance_seed{args.seed}'
    overlay_panels(cases, (args.output_dir or COMPARE / folder).resolve(), args.seed, PANELS, (3.2, 0.4), 3)


def run_one_or_two(args):
    panels = (('NMDAR_World.dat', 'Open NMDAR', 'open_NMDAR'),
              ('NMDAR_Ca_flux_World.dat', r'NMDAR Ca$^{2+}$ entry', 'NMDAR_Ca_entry'),
              ('Ca2_Mushroom_transparent.dat', r'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_Ca'))
    cases = [('1x10_500g', 'One cluster', condition_color('blue', 'dark')),
             ('2x10_100nm_250g', 'Two clusters, 100 nm', condition_color('orange', 'dark'))]
    output = args.output_dir or COMPARE / f'1x10_vs_2x10_seed{args.seed}'
    overlay_panels(cases, output.resolve(), args.seed, panels, (4.2, 0.85), 2)


def run_small_cluster(args):
    output = (args.output_dir or COMPARE / 'small_clusters').resolve()
    families = {1: 'blue', 3: 'purple', 5: 'orange'}
    for n_glu, tone in ((250, 'dark'), (500, 'bright')):
        fig, ax = figure()
        for n_nmdar, family in families.items():
            values = load_trace(trace_path(f'1x{n_nmdar}_{n_glu}g', 1, 'NMDAR_World.dat'))
            ax.plot(values[:, 0] * 1000, values[:, 1], color=condition_color(family, tone),
                    drawstyle='steps-post', label=f'{n_nmdar} NMDAR')
        ax.set(xlabel='Time (ms)', ylabel='Open NMDAR (count)', xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        ax.set_ylim(bottom=-0.05)
        save(fig, output / f'open_NMDAR_{n_glu}g.png', legend_name=None)
        handles = [Line2D([], [], color=condition_color(family, tone), label=f'{count} NMDAR')
                   for count, family in families.items()]
        save_legend(handles, output / f'legend_{n_glu}g.png', figsize=(3.2, 0.4), ncol=3)


# ---- multi-seed batches ----

DISTANCES = (50, 100, 200)
FAMILIES = ('purple', 'orange', 'blue')
GROUPS = (  # receptors per cluster, colour tone
    (10, 'dark'),
    (3, 'bright'),
)


def load_group(batch, seeds, receptor_count):
    expected_time = np.arange(3501) * 1e-5
    traces, warnings = {}, {}
    for distance in DISTANCES:
        case = f'2x{receptor_count}_{distance}nm_250g'
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
            log = work / 'run.log'
            if not log.is_file():
                log = work / 'logs' / f'seed_{seed}.log'
            with log.open('rb') as handle:
                handle.seek(max(0, log.stat().st_size - 65536))
                ending = handle.read().decode(errors='replace')
            if 'Simulation CPU time =' not in ending:
                raise ValueError(f'Missing simulation completion report: {log}')
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
    }
    (output / 'summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for distance in DISTANCES:
        entry = final[str(distance)]['NMDAR_Ca_entry']
        print(f'{receptor_count}/cluster, {distance} nm: final entry '
              f'{entry["mean_final"]:.2f} ± {entry["sample_std_final"]:.2f} (sample SD)')


def common_seeds(batch, receptor_count):
    """Seeds present at every distance of one 2x<receptor_count> series."""
    found = [{int(d.name[5:]) for d in (batch / f'2x{receptor_count}_{distance}nm_250g').glob('seed_[0-9][0-9][0-9][0-9][0-9]')}
             for distance in DISTANCES]
    return sorted(set.intersection(*found))


def run_batch(args):
    batch = args.batch.resolve()
    root = (args.output_dir or COMPARE).resolve()
    # Validate every input before exporting; a group missing any selected seed is skipped.
    loaded = []
    for group in GROUPS:
        seeds = args.seeds or common_seeds(batch, group[0])
        if len(seeds) < 2 or len(set(seeds)) != len(seeds) or any(seed < 1 for seed in seeds):
            print(f'SKIP 2x{group[0]}: needs at least two unique seeds, found {seeds}')
            continue
        try:
            loaded.append((group, seeds, load_group(batch, seeds, group[0])))
        except FileNotFoundError as error:
            print(f'SKIP 2x{group[0]}: {error}')
    for (count, tone), seeds, (time, traces, warnings) in loaded:
        # Default: every seed present; summary.json lists them. Explicit --seeds adds a suffix.
        suffix = '_seeds' + '_'.join(map(str, seeds)) if args.seeds else ''
        plot_group(batch, root / f'2x{count}_distance{suffix}', seeds, count, tone, time, traces, warnings)


# ---- control experiments (run/run_controls.sh), seeds 1-5 ----

CONTROL_PANELS = PANELS + (('RyR_O_World.dat', 'Open RyR (count)', 'open_RyR'),)
CONTROL_SEEDS = (1, 2, 3, 4, 5)
GLUTAMATE = (250, 500, 1000)
DISTANCE_COLORS = dict(zip(DISTANCES, (condition_color(f, 'dark') for f in FAMILIES)))
ONE_CLUSTER_COLORS = ('#8C8C8C', '#555555', '#1A1A1A')
DOSE_STYLES = dict(zip(GLUTAMATE, (':', '--', '-')))


def two_cluster(distance, glu, sites='both'):
    """Case name of a 10+10 condition with glu per site."""
    return f'2x10_{distance}nm_{glu}g' + ('_Aonly' if sites == 'A' else '')


def one_cluster(glu, receptors=10):
    return f'1x{receptors}_{glu}g'


def load_control(runs, label):
    time = np.arange(3501) * 1e-5
    traces = {}
    for seed in CONTROL_SEEDS:
        work = runs / label / f'seed_{seed:05d}'
        for filename, _, panel in CONTROL_PANELS:
            path = work / 'react_data' / f'seed_{seed:05d}' / filename
            if not path.is_file():
                return None, f'{label} seed {seed}'
            values = np.loadtxt(path)
            if values.shape != (3501, 2) or not np.allclose(values[:, 0], time, rtol=0, atol=1e-10):
                return None, f'{label} seed {seed} ({filename} incomplete)'
            traces.setdefault(panel, []).append(values[:, 1])
    return {panel: np.stack(rows) for panel, rows in traces.items()}, None


def control_figures():
    """name -> list of (run label, legend label, color, linestyle)."""
    figures = {'cluster_split': [
        (one_cluster(250), '1x10_250g', ONE_CLUSTER_COLORS[0], '--'),
        (one_cluster(500), '1x10_500g', ONE_CLUSTER_COLORS[1], '--'),
        (one_cluster(500, 20), '1x20_500g', ONE_CLUSTER_COLORS[2], '--'),
    ] + [(two_cluster(d, 250), two_cluster(d, 250), DISTANCE_COLORS[d], '-') for d in DISTANCES]}
    for glu in GLUTAMATE:
        figures[f'2x10_distance_{glu}g'] = [
            (two_cluster(d, glu), f'{d} nm', DISTANCE_COLORS[d], '-') for d in DISTANCES
        ] + [(one_cluster(glu), one_cluster(glu), ONE_CLUSTER_COLORS[1], '--')]
    for d in DISTANCES:
        figures[f'2x10_{d}nm_dose'] = [(two_cluster(d, glu), f'{glu} glu/site', DISTANCE_COLORS[d], DOSE_STYLES[glu])
                                       for glu in GLUTAMATE]
    figures['crosstalk_Aonly'] = [(two_cluster(d, 250, 'A'), f'{d} nm, glu at A only', DISTANCE_COLORS[d], '-')
                                  for d in DISTANCES] + [
        (one_cluster(250), one_cluster(250), ONE_CLUSTER_COLORS[1], '--')]
    return figures


def run_controls(args):
    runs = args.runs.resolve()
    root = (args.output_dir or COMPARE).resolve()
    time_ms = np.arange(3501) * 1e-2
    cache, summary = {}, {}
    for name, cases in control_figures().items():
        missing = []
        for label, *_ in cases:
            if label not in cache:
                cache[label] = load_control(runs, label)
            if cache[label][1]:
                missing.append(cache[label][1])
        if missing:
            print(f'SKIP {name}: missing {", ".join(missing)}')
            continue
        output = root / name
        handles = [Line2D([], [], color=c, linestyle=s, label=l) for _, l, c, s in cases]
        handles.append(Patch(facecolor='gray', alpha=0.16, edgecolor='none',
                             label=f'±1 SD (n={len(CONTROL_SEEDS)})'))
        for _, ylabel, panel in CONTROL_PANELS:
            fig, ax = figure()
            for label, _, color, style in cases:
                samples = cache[label][0][panel]
                mean, std = samples.mean(axis=0), samples.std(axis=0, ddof=1)
                ax.fill_between(time_ms, np.maximum(0, mean - std), mean + std, color=color,
                                alpha=0.16, linewidth=0, step='post')
                ax.plot(time_ms, mean, color=color, linestyle=style, drawstyle='steps-post')
            ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35), ylim=(0, None))
            ax.set_xticks((0, 10, 20, 30))
            save(fig, output / f'{panel}.png', legend_name=None)
        save_legend(handles, output / 'legend.png', figsize=(4.2, 0.25 * (len(handles) + 1)), ncol=1)
        print(f'WROTE {output}')
    # RyR opening is sparse, so each seed is also drawn on its own, in the case's plots/ folder.
    for label, (traces, error) in cache.items():
        if error:
            continue
        fig, ax = figure()
        for row, color in zip(traces['open_RyR'], colors(len(CONTROL_SEEDS))):
            ax.plot(time_ms, row, color=color, drawstyle='steps-post', linewidth=1.2)
        ax.set(xlabel='Time (ms)', ylabel='Open RyR (count)', xlim=(0, 35), ylim=(0, None), title=label)
        ax.set_xticks((0, 10, 20, 30))
        target = (runs / label / 'plots' / 'open_RyR_per_seed.png' if args.output_dir is None
                  else root / 'open_RyR_per_seed' / f'{label}.png')
        save(fig, target, legend_name=None)
        final = lambda panel: [float(v) for v in traces[panel][:, -1]]
        peak = lambda panel: [float(v) for v in traces[panel].max(axis=1)]
        summary[label] = {
            'seeds': list(CONTROL_SEEDS),
            'final_NMDAR_Ca_entry': final('NMDAR_Ca_entry'),
            'peak_open_NMDAR': peak('open_NMDAR'),
            'peak_spine_cytosolic_Ca': peak('spine_cytosolic_Ca'),
            'peak_open_RyR': peak('open_RyR'),
            'seeds_with_RyR_opening': int((traces['open_RyR'].max(axis=1) > 0).sum()),
        }
    root.mkdir(parents=True, exist_ok=True)
    (root / 'controls_summary.json').write_text(json.dumps(dict(
        note='Per-seed values, seeds 1-5. Bands in figures are mean ± sample SD (ddof=1).',
        runs=str(runs),
        conditions=summary), indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('distance')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--receptors', type=int, choices=(3, 10), default=10)
    p.set_defaults(function=run_distance)
    p.add_argument('--output-dir', type=Path)
    p = sub.add_parser('batch')
    p.add_argument('batch', type=Path, nargs='?', default=RUNS, help='Directory containing CASE/seed_NNNNN run folders')
    p.add_argument('--seeds', type=int, nargs='+', help='default: every seed present at all three distances')
    p.add_argument('--output-dir', type=Path)
    p.set_defaults(function=run_batch)
    p = sub.add_parser('controls')
    p.add_argument('--runs', type=Path, default=RUNS, help='CASE/seed_NNNNN run folders')
    p.add_argument('--output-dir', type=Path)
    p.set_defaults(function=run_controls)
    p = sub.add_parser('one-or-two')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--output-dir', type=Path)
    p.set_defaults(function=run_one_or_two)
    p = sub.add_parser('small-cluster')
    p.add_argument('--output-dir', type=Path)
    p.set_defaults(function=run_small_cluster)
    args = parser.parse_args()
    apply_style()
    args.function(args)


if __name__ == '__main__':
    main()
