"""Build an isolated, reproducible run without touching historical models."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parent
CHEMISTRY = ROOT / 'chemistry'
ANALYSIS = ROOT / 'analysis'
DEFAULT_CONDITION = ROOT / 'conditions/2_cluster_100nm.json'
sys.path.insert(0, str(CHEMISTRY))
from calcium_labels import expand_background
from geometry_checks import check_placement


def validate(config):
    assert config['spine_apparatus'] is False
    assert config['time_step_s'] == 1e-7
    n = config['duration_s'] / config['time_step_s']
    assert n >= 1000 and math.isclose(n / 1000, round(n / 1000))
    assert config['viz_every_n_timesteps'] > 0
    assert round(n) % config['viz_every_n_timesteps'] == 0
    assert config['viz_every_n_timesteps'] % 100 == 0
    assert config['voxel_size_um'] > 0
    assert config['channel_protocol']['rest_voltage_mV'] == -67
    assert config['channel_protocol']['magnesium_mM'] == 1
    assert config['calcium_sources'] == {'NMDAR_A': 'Ca_A', 'NMDAR_B': 'Ca_B'}
    assert config['vgcc_populations'] == {'VGCC': 1.0}
    assert list(config['channel_protocol']['nmda']) == ['NMDAR_A', 'NMDAR_B']
    modes = list(config['channel_protocol']['nmda'].values()) + list(config['channel_protocol']['vgcc'].values())
    assert len(set(modes)) == 1 and modes[0] in ('evoked', 'rest')
    assert 0 <= config['waveform_start_s'] < config['duration_s']
    assert len(config['clusters']) == 2
    assert any(c['glutamate'] > 0 for c in config['clusters'])
    seen = set()
    for tag, c in zip('AB', config['clusters']):
        assert c['population'] == 'NMDAR_' + tag
        assert isinstance(c['receptors'], int) and c['receptors'] > 0
        assert isinstance(c['glutamate'], int) and c['glutamate'] >= 0
        assert c['walls'] and not seen.intersection(c['walls'])
        seen.update(c['walls'])
        assert 0 <= c['release_time_s'] < config['duration_s']


def render(config):
    text = expand_background((CHEMISTRY / 'background_no_sa.bngl').read_text())
    kinetics = json.loads((CHEMISTRY / 'nmda_kinetics.json').read_text())
    molecules, rules = [], []
    for population in config['channel_protocol']['nmda']:
        molecules.extend(s.replace('NMDAR', population) for s in kinetics['states'])
        for index, rule in enumerate(kinetics['rules']):
            left, right = rule.replace('NMDAR', population).split('<->')
            product, fwd, rev = re.fullmatch(r'\s*(.*?)\s+(\S+),\s*(\S+)\s*', right).groups()
            rules.extend([f'{population}_gating_{index}_f: {left.strip()} -> {product} {fwd}',
                          f'{population}_gating_{index}_r: {product} -> {left.strip()} {rev}'])
    molecules.extend('VGCC' + suffix for suffix in ('', '_0', '_1', '_2', '_3'))
    text = text.replace('end molecule types', '\n'.join('    ' + s for s in molecules) + '\nend molecule types')
    return text.replace('end reaction rules', '\n'.join('    ' + s for s in rules) + '\nend reaction rules')


def materialize(destination, smoke=False, condition_path=None, glutamate=None):
    config = json.loads((Path(condition_path) if condition_path else DEFAULT_CONDITION).read_text())
    if glutamate is not None:
        for cluster, count in zip(config['clusters'], glutamate):
            cluster['glutamate'] = count
    if smoke:
        config['production_duration_s'] = config['duration_s']
        release = max(c['release_time_s'] for c in config['clusters'])
        config['duration_s'] = round(release + 0.001, 9)
        config['smoke_note'] = 'Baseline plus 1 ms startup/chemistry test; not a complete response.'
    validate(config)
    placement = check_placement(config, ROOT / 'geometry')
    destination = Path(destination).resolve()
    if destination.is_relative_to(ROOT.parents[1] / 'results'):
        raise ValueError('Build inputs under model/calcium_track/snapshots/, not results/.')
    destination.mkdir(parents=True, exist_ok=False)
    for path in (ROOT / 'runtime').glob('*.py'):
        shutil.copy2(path, destination / path.name)
    shutil.copytree(ROOT / 'geometry', destination / 'geometry', ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('geometry_checks.py', 'paths.py', 'batch.py', 'run.py'):
        shutil.copy2(ROOT / name, destination / name)
    for name in ('calcium_labels.py', 'evoked_rates.json', 'resting_rates.json', 'nmda_kinetics.json', 'diffusion.json'):
        shutil.copy2(CHEMISTRY / name, destination / name)
    for path in ANALYSIS.glob('*.py'):
        if path.name != 'compare_spatial.py' and path.name != 'plot.py':
            shutil.copy2(path, destination / path.name)
    (destination / 'condition.json').write_text(json.dumps(config, indent=2) + '\n')
    (destination / 'placement.json').write_text(json.dumps(placement, indent=2) + '\n')
    (destination / 'model.bngl').write_text(render(config))
    hashes = {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(destination.rglob('*')) if p.is_file()}
    (destination / 'input_hashes.json').write_text(json.dumps(hashes, indent=2) + '\n')
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--condition', type=Path, help='Condition JSON from conditions/; default is the 100 nm model.')
    parser.add_argument('--glutamate', type=int, nargs=2, metavar=('A', 'B'),
                        help='Molecules released at each site, e.g. 1000 0 for A only.')
    args = parser.parse_args()
    print(materialize(args.destination, args.smoke, args.condition, args.glutamate))
