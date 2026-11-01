"""Compare seeds without dropping valid zero-response trajectories."""
import sys as _sys
from pathlib import Path as _Path
_here = _Path(__file__).resolve().parent
_sys.path[:0] = [str(_here.parent), str(_here.parent / 'chemistry')]  # source layout; snapshots are flat
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from paths import repository

sys.path.insert(0, str(repository() / 'scripts'))
import plot_style as style


def summarize(batch, seed1):
    records, traces, failures = [], {}, []
    runs = [seed1] + sorted(p for p in batch.glob('seed_*') if p.is_dir())
    for run in runs:
        status = json.loads((run / 'status.json').read_text())
        seed = status['seed']
        if status['state'] != 'complete':
            failures.append(dict(seed=seed, run=str(run), state=status['state']))
            continue
        report = json.loads((run / 'validation.json').read_text())
        data = run / f'react_data/seed_{seed:05d}'
        traces[seed] = {}
        for tag in 'AB':
            free = np.loadtxt(data / f'Ca_{tag}_World.dat')
            opened = np.loadtxt(data / f'NMDAR_{tag}_World.dat')
            entry = np.loadtxt(data / f'NMDAR_{tag}_Ca_flux_World.dat')
            row = dict(seed=seed, source=tag, **report['sources'][tag],
                sampled_conducting_receptor_ms=float(np.sum(opened[:-1,1]*np.diff(opened[:,0]))*1000),
                peak_free_time_ms=float(free[np.argmax(free[:,1]),0]*1000),
                warnings=len(report['warnings']))
            records.append(row)
            traces[seed][tag] = dict(free=free, entry=entry)
    output = batch / 'comparison'
    output.mkdir(exist_ok=False)
    if records:
        with (output / 'seed_metrics.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader(); writer.writerows(records)
    report = dict(seeds_completed=sorted(traces), failed_runs=failures, metrics=records,
        interpretation='Five-seed exploratory comparison; zero-response runs retained. Sampled conducting time is approximate.',
        controls='Same corrected 200 nm condition; only random seed varies. This also resamples initial receptor positions within each patch.')
    (output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    style.apply_style()
    colors = {'A':'#E52B88', 'B':'#169B46'}
    fig, ax = style.figure()
    for tag in 'AB':
        rows = [r for r in records if r['source'] == tag]
        ax.plot([r['seed'] for r in rows], [r['entry_ions'] for r in rows],
                'o-', color=colors[tag], label='Ca_' + tag)
    ax.set(xlabel='Random seed', ylabel='Calcium entry over 35 ms (ions)', ylim=(0,None))
    ax.set_xticks(sorted(traces))
    style.save(fig, output / 'calcium_entry_by_seed'); plt.close(fig)
    seed_colors = ['#234E70','#E69F00','#009E73','#9467BD','#666666']
    for tag in 'AB':
        fig, ax = style.figure()
        for (seed, data), color in zip(sorted(traces.items()), seed_colors):
            values = data[tag]['free']
            ax.plot(values[:,0]*1000, values[:,1], color=color, label=f'Seed {seed}')
        ax.set(xlabel='Time (ms)', ylabel=f'Free Ca_{tag} ions', xlim=(0,35), ylim=(0,None))
        style.save(fig, output / f'free_Ca_{tag}_by_seed'); plt.close(fig)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=Path)
    parser.add_argument('--seed1', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.batch.resolve(), args.seed1.resolve()), indent=2))
