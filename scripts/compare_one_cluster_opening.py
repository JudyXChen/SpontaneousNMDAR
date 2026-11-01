"""Compare the original and calibrated one-cluster seed-1 NMDAR traces."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.plot_style import apply_style, condition_color, figure, save, save_legend  # noqa: E402

ORIGINAL = (ROOT / 'archive/NMDAR_old/results/evoked/conditions/1_cluster_10_NMDAR/'
            'react_data/seed_00001')
CALIBRATED = (ROOT / 'results/calcium_total/evoked/conditions/1_cluster_10_NMDAR/'
              'react_data/seed_00001')
OUTPUT = (ROOT / 'archive/NMDAR_old/comparison/'
          'evoked/compare_1_cluster_10_NMDAR_seed1')


def load_trace(directory: Path, filename: str) -> np.ndarray:
    trace = np.loadtxt(directory / filename)
    if trace.ndim != 2 or trace.shape[1] != 2 or trace[-1, 0] != 0.035:
        raise ValueError(f'Incomplete or malformed trace: {directory / filename}')
    return trace


def sampled_open_runs(open_count: np.ndarray) -> int:
    positive = open_count > 0
    return int(np.count_nonzero(positive & ~np.r_[False, positive[:-1]]))


def summarize(open_trace: np.ndarray, calcium_trace: np.ndarray) -> dict[str, float | int | None]:
    positive_indices = np.flatnonzero(open_trace[:, 1] > 0)
    return {
        'peak_open_NMDAR': int(open_trace[:, 1].max()),
        'first_sampled_open_time_ms': (
            float(open_trace[positive_indices[0], 0] * 1e3)
            if len(positive_indices) else None),
        'last_sampled_open_time_ms': (
            float(open_trace[positive_indices[-1], 0] * 1e3)
            if len(positive_indices) else None),
        'samples_with_open_NMDAR': int(len(positive_indices)),
        'sampled_open_runs': sampled_open_runs(open_trace[:, 1]),
        'integrated_sampled_open_receptor_ms': float(
            np.trapz(open_trace[:, 1], open_trace[:, 0]) * 1e3),
        'cumulative_modeled_NMDAR_Ca_entry_events': int(calcium_trace[-1, 1]),
    }


def main() -> None:
    traces = {}
    metrics = {}
    for name, directory in (('original', ORIGINAL), ('calibrated', CALIBRATED)):
        opened = load_trace(directory, 'NMDAR_World.dat')
        calcium = load_trace(directory, 'NMDAR_Ca_flux_World.dat')
        if not np.array_equal(opened[:, 0], calcium[:, 0]):
            raise ValueError(f'Time-grid mismatch for {name}')
        traces[name] = {'open': opened, 'calcium': calcium}
        metrics[name] = summarize(opened, calcium)
    if not np.array_equal(traces['original']['open'][:, 0],
                          traces['calibrated']['open'][:, 0]):
        raise ValueError('Original and calibrated time grids differ')

    ratios = {
        'integrated_sampled_open_time_ratio_calibrated_over_original': (
            metrics['calibrated']['integrated_sampled_open_receptor_ms']
            / metrics['original']['integrated_sampled_open_receptor_ms']),
        'calcium_entry_event_ratio_calibrated_over_original': (
            metrics['calibrated']['cumulative_modeled_NMDAR_Ca_entry_events']
            / metrics['original']['cumulative_modeled_NMDAR_Ca_entry_events']),
    }
    payload = {
        'comparison': 'one cluster, 10 NMDAR, 500 glutamate, seed 1, 35 ms',
        'sampling_interval_us': 10.0,
        'original': metrics['original'],
        'calibrated': metrics['calibrated'],
        'ratios': ratios,
        'interpretation_limit': (
            'One matched seed is an exploratory realization, not an ensemble '
            'estimate of the kinetic-rate effect.'),
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / 'metrics.json').write_text(json.dumps(payload, indent=2) + '\n')

    time_s = traces['original']['open'][:, 0]
    with (OUTPUT / 'traces.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow([
            'time_s', 'open_original', 'open_calibrated',
            'cumulative_ca_entry_original', 'cumulative_ca_entry_calibrated'])
        for index, time in enumerate(time_s):
            writer.writerow([
                f'{time:.8g}',
                f'{traces["original"]["open"][index, 1]:.8g}',
                f'{traces["calibrated"]["open"][index, 1]:.8g}',
                f'{traces["original"]["calcium"][index, 1]:.8g}',
                f'{traces["calibrated"]["calcium"][index, 1]:.8g}',
            ])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    apply_style()
    trace_colors = {
        'original': condition_color('blue', 'bright'),
        'calibrated': condition_color('blue', 'dark'),
    }
    linestyles = {'original': 'solid', 'calibrated': 'solid'}
    linewidths = {'original': 1.6, 'calibrated': 1.6}
    plot_specs = (
        ('open_NMDAR', 'open', 'Open NMDAR'),
        ('cumulative_NMDAR_Ca_entry', 'calcium',
         r'NMDAR Ca$^{2+}$ entry'),
    )
    for filename, trace_name, ylabel in plot_specs:
        fig, ax = figure()
        for name in ('calibrated', 'original'):
            trace = traces[name][trace_name]
            ax.plot(trace[:, 0] * 1e3, trace[:, 1], color=trace_colors[name],
                    linestyle=linestyles[name], linewidth=linewidths[name],
                    drawstyle='steps-post', label=name.capitalize())
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if trace_name == 'open':
            ax.set_ylim(-0.05, 1.05)
            ax.set_yticks((0, 1))
        save(fig, OUTPUT / f'{filename}.png', legend_name=None)

    handles = [Line2D([0], [0], color=trace_colors[name], linestyle=linestyles[name],
                      linewidth=linewidths[name],
                      label=name.capitalize())
               for name in ('original', 'calibrated')]
    save_legend(handles, OUTPUT / 'legend.png', figsize=(3, 0.65), ncol=2)
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
