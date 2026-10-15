"""Validate complete traces, conserved populations, and additive observables."""
import json
import math
from pathlib import Path


def validate_run(work, seed, config):
    work = Path(work)
    data = work / 'react_data' / f'seed_{seed:05d}'
    iterations = round(config['duration_s'] / config['time_step_s'])
    cache = {}

    def trace(name):
        if name not in cache:
            rows = [list(map(float, line.split())) for line in (data / (name + '.dat')).read_text().splitlines()]
            if len(rows) != iterations // 100 + 1:
                raise ValueError(f'{name}: incomplete time grid ({len(rows)} samples)')
            for i, row in enumerate(rows):
                if len(row) != 2 or not all(math.isfinite(v) for v in row):
                    raise ValueError(f'{name}: invalid numeric sample')
                if not math.isclose(row[0], i * 100 * config['time_step_s'], abs_tol=1e-10, rel_tol=0):
                    raise ValueError(f'{name}: mismatched sampling time')
            cache[name] = [row[1] for row in rows]
        return cache[name]

    metrics = {'condition': config['condition'], 'seed': seed, 'duration_s': config['duration_s']}
    release_time = min(c['release_time_s'] for c in config['clusters'])
    release_index = round(release_time / (100 * config['time_step_s']))
    for kind, populations in [('NMDAR', config['channel_protocol']['nmda']),
                               ('VGCC', config['channel_protocol']['vgcc'])]:
        for suffix in ('_World', '_total_World', '_Ca_flux_World'):
            total = trace(kind + suffix)
            parts = [trace(population + suffix) for population in populations]
            if any(not math.isclose(a, sum(b), abs_tol=1e-8, rel_tol=0) for a, b in zip(total, zip(*parts))):
                raise ValueError(kind + suffix + ': combined count does not equal population sum')
        for population in populations:
            counts = trace(population + '_total_World')
            expected = (sum(c['receptors'] for c in config['clusters'] if c['population'] == population)
                        if kind == 'NMDAR' else counts[0])
            if any(v != expected for v in counts):
                raise ValueError(population + ': receptor total was not conserved')
            opened = trace(population + '_World')
            if any(v < 0 or v > expected for v in opened):
                raise ValueError(population + ': open count exceeds receptor total')
            entry = trace(population + '_Ca_flux_World')
            if any(b < a for a, b in zip(entry, entry[1:])):
                raise ValueError(population + ': calcium-entry count decreased')
            step = 100 * config['time_step_s']
            integral = sum((a+b)*.5*step for a,b in zip(opened[release_index:], opened[release_index+1:]))
            metrics[population] = dict(total_receptors=expected, peak_open_post_release=max(opened[release_index:]),
                                       open_receptor_seconds_post_release=integral,
                                       calcium_entry_firings_post_release=entry[-1]-entry[release_index])
    for alias, canonical in [('VSCC_World','VGCC_World'),('VSCC_Ca_flux_World','VGCC_Ca_flux_World')]:
        if trace(alias) != trace(canonical):
            raise ValueError('Legacy VGCC alias differs: ' + alias)
    trace('Ca2_Mushroom_transparent')
    log = (work / 'run.log').read_text()
    for i, c in enumerate(config['clusters'], 1):
        expected = f'Released {c["receptors"]} {c["population"]}_0 from "NMDAR_cluster_{i}" at iteration 0.'
        if expected not in log:
            raise ValueError('Missing expected receptor release: ' + expected)
        iteration = round(c['release_time_s'] / config['time_step_s'])
        expected = f'Released {c["glutamate"]} glu from "glutamate_cluster_{i}" at iteration {iteration}.'
        if expected not in log:
            raise ValueError('Missing expected glutamate release: ' + expected)
    metrics['warnings'] = [line for line in log.splitlines() if line.lower().startswith('warning:')]
    (work / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    return metrics
