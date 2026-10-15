"""Diagnostic only: initially open channels exercise nonzero flux arithmetic.

This builds an isolated snapshot and changes only that snapshot's initial
channel states. It is not a production condition or a physiological test.
Usage: MCELL_PATH=... PYTHON=... "$PYTHON" -B check_nonzero_flux.py NEW_OUTPUT
"""
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build import load_condition, materialize
from batch import suite_leaves


def main():
    out = Path(sys.argv[1]).resolve()
    config = load_condition(suite_leaves('mixed_distance')[0])
    config['duration_s'] = 1e-5
    config['waveform_start_s'] = 0.
    for cluster in config['clusters']:
        cluster['release_time_s'] = 0.
    config['diagnostic_only'] = 'All NMDAR and VGCC channels initially open; not a production condition.'
    materialize(config, out)
    path = out / 'instantiation.py'
    code = path.read_text()
    assert code.count("population + '_0'") == 1
    assert code.count("cluster['population'] + '_0'") == 1
    code = code.replace("population + '_0'", 'population')
    code = code.replace("cluster['population'] + '_0'", "cluster['population']")
    path.write_text(code)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    for key in list(env):
        if key.startswith('MCELL_') and key != 'MCELL_PATH':
            del env[key]
    with (out / 'run.log').open('w') as log:
        subprocess.run([sys.executable, '-B', 'model.py', '-seed', '1'], cwd=out,
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    data = out / 'react_data/seed_00001'

    def values(name):
        rows = [list(map(float, line.split())) for line in (data / (name + '.dat')).read_text().splitlines()]
        assert len(rows) == 2 and rows[-1][0] == 1e-5, name
        return [r[1] for r in rows]

    report = {}
    for population in ('eNMDAR', 'sNMDAR', 'eVGCC'):
        total = values(population + '_total_World')
        assert total[0] == total[1] and total[0] > 0
        entry = values(population + '_Ca_flux_World')
        assert entry[-1] > 0, population
        report[population] = {'total': total[-1], 'entry_firings': entry[-1]}
    for suffix in ('_World', '_total_World', '_Ca_flux_World'):
        assert values('NMDAR' + suffix) == [a+b for a,b in zip(values('eNMDAR' + suffix), values('sNMDAR' + suffix))]
        assert values('VGCC' + suffix) == values('eVGCC' + suffix)
    assert values('VSCC_Ca_flux_World') == values('eVGCC_Ca_flux_World')
    (out / 'diagnostic_result.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: nonzero population fluxes, aggregate sums, and conserved channel totals.')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
