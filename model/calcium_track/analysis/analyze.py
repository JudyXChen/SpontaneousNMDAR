"""Validate a completed run and bin saved calcium positions into virtual voxels.

Voxels are analysis bins, never MCell objects or diffusion barriers. Rows contain
particle counts, not concentrations (boundary voxels have partial cytosolic volume).
Both free and bound calcium are exported, with carrier identity retained.
"""
import sys as _sys
from pathlib import Path as _Path
_here = _Path(__file__).resolve().parent
_sys.path[:0] = [str(_here.parent), str(_here.parent / 'chemistry')]  # source layout; snapshots are flat
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
from paths import inputs_for, verify_inputs, repository
import re
import numpy as np
from calcium_labels import CARRIERS, labeled


def analyze(run, seed=1):
    config = json.loads((inputs_for(run) / 'condition.json').read_text())
    dt = config['time_step_s']
    iterations = round(config['duration_s'] / dt)
    folder = run / f'react_data/seed_{seed:05d}'
    viz = run / f'viz_data/seed_{seed:05d}'
    verify_inputs(inputs_for(run))
    traces = {}
    expected_times = np.arange(iterations // 100 + 1) * 100 * dt
    for path in folder.glob('*.dat'):
        data = np.loadtxt(path, ndmin=2)
        assert data.shape == (len(expected_times), 2), str(path)
        assert np.allclose(data[:, 0], expected_times, rtol=0, atol=1e-10), str(path)
        assert np.all(data[:, 1] >= 0), str(path)
        traces[path.stem] = data[:, 1]

    def get(name):
        return traces[name + '_World']

    for tag, cluster in zip('AB', config['clusters']):
        assert np.all(get('NMDAR_' + tag + '_total') == cluster['receptors'])
        assert np.all(get('NMDAR_' + tag) <= cluster['receptors'])
        flux = get('NMDAR_' + tag + '_Ca_flux')
        assert np.all(np.diff(flux) >= 0)
        sinks = sum(get(s + '_' + tag) for s in ('Ca2Null', 'rxn_PMCA', 'rxn_NCX'))
        inventory = get('Ca_' + tag + '_inventory')
        assert inventory[0] == 0
        # Any source conversion or uncounted particle loss fails this identity.
        assert np.array_equal(flux, inventory + sinks), f'Ca_{tag} source balance'
    assert np.array_equal(get('NMDAR_Ca_flux'), get('NMDAR_A_Ca_flux') + get('NMDAR_B_Ca_flux'))
    assert np.array_equal(get('Ca_free_total'), get('Ca2') + get('Ca_A') + get('Ca_B'))
    release_sample = round(min(c['release_time_s'] for c in config['clusters'] if c['glutamate'] > 0) / dt) // 100
    assert get('glu')[release_sample] == sum(c['glutamate'] for c in config['clusters'])
    assert np.all(get('glu')[:release_sample] == 0), 'glutamate present before release'
    for name in ('Bm', 'Bf', 'PMCA', 'NCX', 'PLC', 'VGCC'):
        values = get(name + '_total')
        assert np.all(values == values[0]), name + ' capacity conservation'

    data_model = json.loads(next(viz.glob('Scene.data_model.*.json')).read_text())['mcell']
    assert not any(o['name'] == 'SA' for o in data_model['geometrical_objects']['object_list'])
    colors = {m['mol_name']: m['display']['color'] for m in data_model['define_molecules']['molecule_list']}
    assert colors['Ca_A'] == [1.0, 0.12, 0.55] and colors['Ca_B'] == [0.1, 0.9, 0.2]
    files = sorted(viz.glob('Scene.ascii.*.dat'))
    expected_iterations = list(range(0, iterations + 1, config['viz_every_n_timesteps']))
    assert [int(p.name.split('.')[-2]) for p in files] == expected_iterations
    species = {c: ('background', 'free' if c == 'Ca2' else 'bound') for c in CARRIERS}
    species.update({labeled(c, tag): (tag, 'free' if c == 'Ca2' else 'bound')
                    for tag in 'AB' for c in CARRIERS})
    size = config['voxel_size_um']
    output = run / 'spatial'
    output.mkdir(exist_ok=True)
    # World coordinates preserve the inherited spine and dendrite geometry.
    # Origin is fixed at (0,0,0); negative indices are valid.
    with gzip.open(output / 'calcium_voxels.csv.gz', 'wt', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['iteration','time_s','ix','iy','iz','source','state','carrier','count'])
        for path, iteration in zip(files, expected_iterations):
            bins, counts = Counter(), Counter()
            for line in path.read_text().splitlines():
                fields = line.split()
                name = fields[0]
                if name not in species:
                    continue
                xyz = np.array([float(s) for s in fields[2:5]])
                index = tuple(np.floor(xyz / size).astype(int))
                bins[(*index, *species[name], name)] += 1
                counts[name] += 1
            for name in species:
                assert counts[name] == get(name)[iteration // 100], f'Viz/count mismatch: {path.name} {name}'
            for key, count in sorted(bins.items()):
                writer.writerow([iteration, f'{iteration * dt:.9g}', *key, count])
    metadata = dict(voxel_size_um=size, origin_um=[0,0,0], frame_interval_s=config['viz_every_n_timesteps']*dt,
                    frames=len(files), coordinates='Inherited geometry world coordinates, micrometers.',
                    extent='All recorded calcium in spine and attached dendrite; keep spatial coordinates to select an ROI.',
                    missing_rows='Zero occupancy at the listed saved frame; frame list is explicit below.',
                    frame_iterations=expected_iterations,
                    quantities='Particle counts, not concentration; no correction for partial voxel cytosolic volume.',
                    labels='A/B identify the NMDAR entry cluster, not the glutamate release site. Bound carriers hold one calcium each.',
                    temporal_limit='Only saved instants are recoverable; no claim about paths between frames.')
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    log = (run / 'simulation.log').read_text()
    warnings = [line for line in log.splitlines() if re.search('warning|missed|probability.*(?:exceed|greater)', line, re.I)]
    metrics = dict(validation='passed', duration_s=config['duration_s'], samples=len(expected_times),
                   frames=len(files), source_balance='entry = free + bound + explicit sinks at every count sample',
                   warnings=warnings,
                   sources={tag: dict(entry_ions=int(get('NMDAR_' + tag + '_Ca_flux')[-1]),
                       final_free=int(get('Ca_' + tag)[-1]), final_inventory=int(get('Ca_' + tag + '_inventory')[-1]),
                       peak_free=int(get('Ca_' + tag).max())) for tag in 'AB'})
    (run / 'validation.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps(metrics, indent=2))
    return metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()
    analyze(args.run_directory.resolve(), args.seed)
