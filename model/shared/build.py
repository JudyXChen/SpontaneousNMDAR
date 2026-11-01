"""Resolve small condition files into self-contained MCell input snapshots."""
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

COMMON = Path(__file__).resolve().parent
MODEL_ROOT = COMMON.parent


def read_json(path):
    return json.loads(path.read_text())


def load_condition(leaf, smoke=False):
    leaf = Path(leaf).resolve()
    config = dict(read_json(leaf / 'instantiation.json'), **read_json(leaf / 'parameters.json'))
    config['channel_protocol'] = read_json(leaf.parent / 'channels.json')
    config['condition'] = str(leaf.relative_to(MODEL_ROOT))
    validate(config)
    if smoke:
        config['production_timeline'] = {k: copy.deepcopy(config[k]) for k in
                                         ('clusters', 'duration_s', 'waveform_start_s')}
        config['duration_s'] = 0.0002
        config['waveform_start_s'] = 0.0001
        for cluster in config['clusters']:
            cluster['release_time_s'] = 0.0001
        config['smoke_note'] = 'Compressed 0.1 ms baseline + 0.1 ms response; not production baseline validation.'
    validate(config)
    return config


def validate(config):
    dt = config['time_step_s']
    if dt != 1e-7:
        raise ValueError('These protocols require time_step_s = 1e-7')
    n = config['duration_s'] / dt
    if n < 100 or not math.isclose(n / 100, round(n / 100), abs_tol=1e-8):
        raise ValueError('Duration must be a positive multiple of 100 timesteps')
    protocol = config['channel_protocol']
    if protocol['rest_voltage_mV'] != -67 or protocol['magnesium_mM'] != 1:
        raise ValueError('Resting rate data currently encode -67 mV and 1 mM Mg only')
    if not protocol['nmda'] or not protocol['vgcc']:
        raise ValueError('A channel population must be declared for each channel class')
    for population, mode in {**protocol['nmda'], **protocol['vgcc']}.items():
        expected = {'sNMDAR': 'rest', 'eNMDAR': 'evoked', 'sVGCC': 'rest', 'eVGCC': 'evoked'}
        if population not in expected or mode != expected[population]:
            raise ValueError('Population identity must agree with its drive: ' + population)
    if set(config['vgcc_populations']) != set(protocol['vgcc']):
        raise ValueError('VGCC placement and channel definitions must agree')
    fractions = list(config['vgcc_populations'].values())
    if any(x <= 0 for x in fractions) or not math.isclose(sum(fractions), 1):
        raise ValueError('VGCC fractions must sum to one; do not duplicate the background density')
    seen = set()
    for c in config['clusters']:
        if c['population'] not in protocol['nmda']:
            raise ValueError('Unknown NMDAR population')
        for key in ('receptors', 'glutamate'):
            if type(c[key]) is not int or c[key] < 1:
                raise ValueError(key + ' must be a positive integer')
        if not c['walls'] or any(type(w) is not int or w < 0 for w in c['walls']):
            raise ValueError('Cluster walls must be nonempty lists of nonnegative integers')
        if seen.intersection(c['walls']) or len(set(c['walls'])) != len(c['walls']):
            raise ValueError('Cluster wall regions must not overlap or contain duplicates')
        seen.update(c['walls'])
        if not 0 <= c['release_time_s'] < config['duration_s']:
            raise ValueError('Glutamate must be released within the simulation')
        if not math.isclose(c['release_time_s']/dt, round(c['release_time_s']/dt), abs_tol=1e-6):
            raise ValueError('Release times must align with the timestep')
    if not 0 <= config['waveform_start_s'] < config['duration_s']:
        raise ValueError('Waveform start must be within the simulation')


def render_bngl(config):
    text = (COMMON / 'background.bngl').read_text()
    kinetics = read_json(COMMON / 'nmda_kinetics.json')
    molecules, rules = [], []
    for population in config['channel_protocol']['nmda']:
        molecules.extend(name.replace('NMDAR', population) for name in kinetics['states'])
        for i, rule in enumerate(kinetics['rules'], 1):
            renamed = re.sub(r'\bNMDAR(?:_[A-Za-z0-9]+)?\b',
                             lambda m: m[0].replace('NMDAR', population), rule)
            # MCell assigns the same label to both directions of a named
            # reversible BNGL rule. Give the two directions unique names.
            left, right = renamed.split('<->')
            match = re.fullmatch(r'\s*(.*?)\s+(\S+),\s*(\S+)\s*', right)
            product, forward, reverse = match.groups()
            rules.append(f'{population}_gating_{i}_f: {left.strip()} -> {product} {forward}')
            rules.append(f'{population}_gating_{i}_r: {product} -> {left.strip()} {reverse}')
    for population in config['channel_protocol']['vgcc']:
        molecules.extend(population + suffix for suffix in ('', '_0', '_1', '_2', '_3'))
    text = text.replace('end molecule types', '\n'.join('    ' + x for x in molecules) + '\nend molecule types')
    text = text.replace('end reaction rules', '\n'.join('    ' + x for x in rules) + '\nend reaction rules')
    text = re.sub(r'(?m)^(\s*n_NMDAR)\s+\S+',
                  lambda m: m[1] + ' ' + str(sum(c['receptors'] for c in config['clusters'])), text)
    text = re.sub(r'(?m)^(\s*n_glu)\s+\S+',
                  lambda m: m[1] + ' ' + str(config['clusters'][0]['glutamate']), text)
    return text


def materialize(config, destination):
    """Never overwrite a snapshot or copy old simulation outputs/checkpoints."""
    validate(config)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    for path in (COMMON / 'runtime').glob('*.py'):
        shutil.copy2(path, destination / path.name)
    shutil.copytree(COMMON / 'geometry', destination / 'geometry',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('evoked_rates.json', 'resting_rates.json', 'nmda_kinetics.json'):
        shutil.copy2(COMMON / name, destination / name)
    (destination / 'condition.json').write_text(json.dumps(config, indent=2) + '\n')
    (destination / 'model.bngl').write_text(render_bngl(config))
    return {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(destination.rglob('*')) if p.is_file()}
