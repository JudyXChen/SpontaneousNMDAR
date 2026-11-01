"""Regenerate calcium_total figures from recorded data. No simulations are run.

Run with the Blender Python or another Python with NumPy and Matplotlib.
Existing figures are overwritten in place; nothing else is removed.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results/calcium_total'
MODELS = ROOT / 'model/calcium_total'


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    os.environ.setdefault('MPLBACKEND', 'Agg')
    tasks = []
    if (RESULTS / 'conductance/comparison.csv').exists():
        tasks.append(('plot_conductance.py', '--from-csv', '--output-dir', RESULTS / 'conductance'))
    for protocol in ('evoked', 'spont'):
        for output in sorted((RESULTS / protocol / 'conditions').iterdir()):
            model = MODELS / protocol / output.name
            if (output / 'ca_er.png').exists():
                tasks.append(('plot_condition.py', 'calcium', model))
            if (output / 'flux_nmdar.png').exists():
                tasks.append(('plot_condition.py', 'flux', model))
            if protocol == 'evoked' and (output / 'psd_layout.png').exists():
                tasks.append(('plot_condition.py', 'traces', model))
    tasks += [('plot_evoked_compare.py', 'one-or-two'),
              ('plot_evoked_compare.py', 'distance'),
              ('plot_evoked_compare.py', 'distance', '--receptors', '3'),
              ('plot_evoked_compare.py', 'small-cluster'),
              ('plot_spont_compare.py', '--raw-only', '--output-dir', RESULTS / 'spont/compare_1cluster')]
    for index, (script, *arguments) in enumerate(tasks, 1):
        print(f'[{index}/{len(tasks)}] {script} {" ".join(map(str, arguments))}', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'scripts' / script), *map(str, arguments)], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
