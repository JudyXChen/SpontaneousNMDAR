"""Run a prepared snapshot, preserve status/logs, then validate and voxelize."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from paths import inputs_for, verify_inputs

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / 'analysis' if (HERE / 'analysis').is_dir() else HERE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--inputs', type=Path, required=True)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    inputs = args.inputs.resolve()
    verify_inputs(inputs)
    run.mkdir(parents=True, exist_ok=False)
    (run / 'input_reference.json').write_text(json.dumps({
        'input_directory': os.path.relpath(inputs, run),
        'seed': args.seed, 'layout': 'Inputs outside results; shared immutable snapshot.'
    }, indent=2) + '\n')
    if (run / 'status.json').exists():
        raise FileExistsError('Run already started; prepare a new snapshot instead.')
    started = time.monotonic()
    status = dict(state='running', seed=args.seed, start_utc=datetime.now(timezone.utc).isoformat(),
                  python=sys.executable, mcell_path=os.environ.get('MCELL_PATH'), pid=os.getpid())

    def save():
        status['elapsed_s'] = time.monotonic() - started
        temporary = run / 'status.tmp'
        temporary.write_text(json.dumps(status, indent=2) + '\n')
        temporary.replace(run / 'status.json')

    save()
    try:
        with (run / 'simulation.log').open('w') as log:
            process = subprocess.run([sys.executable, '-B', str(inputs / 'model.py'), '-seed', str(args.seed)],
                                     cwd=run, stdout=log, stderr=subprocess.STDOUT)
        status['simulation_returncode'] = process.returncode
        status['simulation_elapsed_s'] = time.monotonic() - started
        if process.returncode:
            raise RuntimeError('MCell failed; see simulation.log')
        status['state'] = 'analyzing'
        save()
        with (run / 'analysis.log').open('w') as log:
            subprocess.run([sys.executable, '-B', str(ANALYSIS / 'analyze.py'),
                            str(run), '--seed', str(args.seed)], stdout=log,
                           stderr=subprocess.STDOUT, check=True)
            subprocess.run([sys.executable, '-B', str(ANALYSIS / 'audit.py'),
                            str(run), '--output', str(run / 'audit.json')], stdout=log,
                           stderr=subprocess.STDOUT, check=True)
            subprocess.run([sys.executable, '-B', str(ANALYSIS / 'export_clean_viz.py'),
                            str(run), '--seed', str(args.seed)], stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        status['state'] = 'complete'
    except BaseException as error:
        status['state'] = 'failed'
        status['error'] = repr(error)
        raise
    finally:
        status['end_utc'] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
