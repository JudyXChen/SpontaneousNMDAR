"""Plot calcium_track count traces as seed means with sample SD bands.

  plot_calcium_traces.py OUTPUT BATCH                 per-cluster NMDAR and calcium panels
  plot_calcium_traces.py OUTPUT NAME=BATCH NAME=BATCH conditions overlaid, clusters summed

Each BATCH holds seed_NNNNN run folders. Free calcium is counted in the spine (head,
neck and root), and "total" adds Ca_A, Ca_B and unlabeled Ca2. Channels flicker faster
than the panel can show, so open counts are smoothed over 1 ms; traces.csv keeps raw means.
"""
import argparse
import csv
from pathlib import Path
import sys

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as style

PINK, GREEN = '#E52B88', '#169B46'
OPEN_LABEL, ENTRY_LABEL = 'Open NMDAR (count, 1 ms mean)', 'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)'
TOTAL_LABEL = 'Spine cytosolic free Ca$^{2+}$ (ions)'


def load_batch(batch):
    """Return time and {trace name: array of shape (seeds, samples)}."""
    runs = sorted(path for path in batch.glob('seed_*') if path.is_dir())
    if not runs:
        raise ValueError(f'No seed folders in {batch}')
    series, time = {}, None
    for run in runs:
        data = next(path for path in (run / 'react_data').iterdir() if path.is_dir())
        for name, file_name in (('open_A', 'NMDAR_A_World'), ('open_B', 'NMDAR_B_World'),
                                ('entry_A', 'NMDAR_A_Ca_flux_World'), ('entry_B', 'NMDAR_B_Ca_flux_World'),
                                ('free_A', 'Ca_A_Mushroom_transparent'), ('free_B', 'Ca_B_Mushroom_transparent'),
                                ('free_background', 'Ca2_Mushroom_transparent')):
            values = np.loadtxt(data / (file_name + '.dat'))
            if time is None:
                time = values[:, 0]
            elif not np.array_equal(time, values[:, 0]):
                raise ValueError(f'Time grid differs: {data / file_name}')
            series.setdefault(name, []).append(values[:, 1])
    series = {name: np.stack(values) for name, values in series.items()}
    series['free_total'] = series['free_A'] + series['free_B'] + series['free_background']
    series['open_all'] = series['open_A'] + series['open_B']
    series['entry_all'] = series['entry_A'] + series['entry_B']
    return time * 1000, series, len(runs)


def smooth(samples, width=100):
    """Moving mean per seed over width samples (10 us each), shorter at the ends."""
    kernel = np.ones(width)
    weights = np.convolve(np.ones(samples.shape[1]), kernel, mode='same')
    return np.stack([np.convolve(row, kernel, mode='same') / weights for row in samples])


def panel(time, curves, ylabel, path, band_note):
    fig, ax = style.figure()
    handles = []
    for samples, label, color in curves:
        mean, std = samples.mean(axis=0), samples.std(axis=0, ddof=1)
        ax.fill_between(time, np.maximum(0, mean - std), mean + std, color=color, alpha=0.16,
                        linewidth=0, step='post')
        ax.plot(time, mean, color=color, drawstyle='steps-post', label=label)
        handles.append(Line2D([], [], color=color, label=label))
    handles.append(Patch(facecolor='gray', alpha=0.16, edgecolor='none', label=band_note))
    ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, time[-1]), ylim=(0, None))
    ax.set_xticks(np.arange(0, time[-1] + 1e-9, 10))
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    style.save(fig, path, legend_name=None)
    return handles


def write_csv(path, time, columns):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_ms'] + [f'{name}_{kind}' for name in columns for kind in ('mean', 'std')])
        rows = [np.column_stack((values.mean(axis=0), values.std(axis=0, ddof=1))) for values in columns.values()]
        for index, moment in enumerate(time):
            writer.writerow([f'{moment:.6g}'] + [f'{value:.6g}' for row in rows for value in row[index]])
    print(f'Saved -> {path}')


def one_batch(output, batch):
    time, s, n = load_batch(batch)
    note = f'±1 SD (n={n})'
    clusters = lambda kind: [(s[f'{kind}_A'], 'Cluster A', PINK), (s[f'{kind}_B'], 'Cluster B', GREEN)]
    panel(time, [(smooth(values), label, color) for values, label, color in clusters('open')], OPEN_LABEL, output / 'open_NMDAR.png', note)
    handles = panel(time, clusters('entry'), ENTRY_LABEL, output / 'NMDAR_Ca_entry.png', note)
    style.save_legend(handles, output / 'legend.png', figsize=(3.2, 0.65), ncol=2)
    sources = [(s['free_A'], 'Ca_A', PINK), (s['free_B'], 'Ca_B', GREEN),
               (s['free_background'], 'Unlabeled Ca2', style.shade('blue', 'dark'))]
    handles = panel(time, sources, 'Spine free Ca$^{2+}$ by source (ions)', output / 'spine_free_Ca_by_source.png', note)
    style.save_legend(handles, output / 'sources_legend.png', figsize=(3.2, 0.65), ncol=2)
    panel(time, [(s['free_total'], 'Total', 'black')], TOTAL_LABEL, output / 'spine_cytosolic_Ca.png', note)
    write_csv(output / 'traces.csv', time, s)


def compare(output, pairs):
    loaded = [(name, *load_batch(batch)) for name, batch in pairs]
    time = loaded[0][1]
    if any(not np.array_equal(time, item[1]) for item in loaded):
        raise ValueError('Conditions use different time grids')
    colors = style.colors(len(loaded))
    note = '±1 SD (n=' + ', '.join(str(item[3]) for item in loaded) + ')'
    for key, ylabel, stem in (('open_all', OPEN_LABEL, 'open_NMDAR'), ('entry_all', ENTRY_LABEL, 'NMDAR_Ca_entry'),
                              ('free_total', TOTAL_LABEL, 'spine_cytosolic_Ca')):
        curve = smooth if key == 'open_all' else (lambda values: values)
        handles = panel(time, [(curve(item[2][key]), item[0], color) for item, color in zip(loaded, colors)],
                        ylabel, output / f'{stem}.png', note)
    style.save_legend(handles, output / 'legend.png', figsize=(3.2, 0.65), ncol=2)
    write_csv(output / 'traces.csv', time, {f'{key}_{item[0]}': item[2][key] for item in loaded
                                            for key in ('open_all', 'entry_all', 'free_total')})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output', type=Path)
    parser.add_argument('batches', nargs='+', help='BATCH, or NAME=BATCH for each compared condition')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style.apply_style()
    if len(args.batches) == 1 and '=' not in args.batches[0]:
        one_batch(args.output, Path(args.batches[0]))
    else:
        compare(args.output, [(item.split('=', 1)[0], Path(item.split('=', 1)[1])) for item in args.batches])


if __name__ == '__main__':
    main()
