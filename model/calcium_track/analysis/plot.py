"""PNG previews from a completed source-calcium run; no simulation changes."""
import sys as _sys
from pathlib import Path as _Path
_here = _Path(__file__).resolve().parent
_sys.path[:0] = [str(_here.parent), str(_here.parent / 'chemistry')]  # source layout; snapshots are flat
import argparse
import ast
import json
from pathlib import Path
from paths import inputs_for, verify_inputs, repository
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

sys.path.insert(0, str(repository() / 'scripts'))
import plot_style as style


def plot(run, seed):
    assert json.loads((run / 'validation.json').read_text())['validation'] == 'passed'
    style.apply_style()
    data = run / f'react_data/seed_{seed:05d}'
    output = run / 'plots'
    palette = {'A': '#E52B88', 'B': '#169B46'}
    traces = {tag: np.loadtxt(data / f'Ca_{tag}_World.dat') for tag in 'AB'}
    for kind, pattern, ylabel in [('free_calcium', 'Ca_{}_World.dat', 'Free calcium ions'),
                                  ('calcium_entry', 'NMDAR_{}_Ca_flux_World.dat', 'Cumulative calcium entry (ions)')]:
        fig, ax = style.figure()
        for tag in 'AB':
            values = np.loadtxt(data / pattern.format(tag))
            ax.plot(values[:, 0]*1000, values[:, 1], color=palette[tag], label='Cluster ' + tag)
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0,traces['A'][-1,0]*1000), ylim=(0,None))
        style.save(fig, output / kind)
        plt.close(fig)
    config = json.loads((inputs_for(run) / 'condition.json').read_text())
    step = config['viz_every_n_timesteps'] // 100
    joint = traces['A'][::step, 1] + traces['B'][::step, 1]
    frame = int(np.argmax(joint))
    iteration = frame * config['viz_every_n_timesteps']
    time_ms = iteration * config['time_step_s'] * 1000
    path = next(p for p in (run / f'viz_data/seed_{seed:05d}').glob('Scene.ascii.*.dat')
                if int(p.name.split('.')[-2]) == iteration)
    points = {tag: [] for tag in 'AB'}
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields[0] in ('Ca_A','Ca_B'):
            points[fields[0][-1]].append([float(s) for s in fields[2:5]])
    tree = ast.parse((inputs_for(run) / 'geometry/Mushroom_transparent.py').read_text())
    constants = {n.targets[0].id: ast.literal_eval(n.value) for n in tree.body
                 if isinstance(n, ast.Assign) and isinstance(n.value, ast.List)}
    vertices = np.asarray(constants['Mushroom_transparent_vertex_list'])
    walls = np.asarray(constants['Mushroom_transparent_wall_list'])
    edges = np.concatenate([walls[:, [0,1]], walls[:, [1,2]], walls[:, [2,0]]])
    fig, ax = style.figure()
    ax.add_collection(LineCollection(vertices[edges][:, :, [0,2]], color='#999999', linewidths=0.2, alpha=0.15))
    for tag in 'AB':
        xyz = np.asarray(points[tag]).reshape(-1,3)
        ax.scatter(xyz[:,0], xyz[:,2], s=8, color=palette[tag], alpha=0.7, linewidths=0, label='Ca_' + tag)
    ax.set(xlabel='x (µm)', ylabel='z (µm)', xlim=(-0.7,0.7), ylim=(0.35,1.75))
    ax.set_aspect('equal')
    style.save(fig, output / f'free_calcium_xz_{time_ms:06.3f}ms.png')
    plt.close(fig)
    (output / 'preview_metadata.json').write_text(json.dumps(dict(
        time_ms=time_ms, iteration=iteration,
        selection='Saved frame with maximum combined free Ca_A + Ca_B across the full cytosolic domain.',
        projection='x-z projection; y collapsed. View is cropped to the spine; no bound calcium shown.',
        counts_world={tag: len(points[tag]) for tag in 'AB'}), indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()
    plot(args.run_directory.resolve(), args.seed)
