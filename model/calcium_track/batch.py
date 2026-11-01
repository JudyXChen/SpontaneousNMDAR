"""Run a batch plan's seeds sequentially and summarize recorded outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
from paths import repository, verify_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--workers', type=int, help='Seeds to run at once; default comes from the plan.')
    args = parser.parse_args()
    repo = repository()
    plan = json.loads(args.plan.read_text())
    inputs, output = repo / plan['inputs'], repo / plan['outputs']
    verify_inputs(inputs)
    output.mkdir(parents=True, exist_ok=False)
    (output / 'input_reference.json').write_text(json.dumps(dict(
        input_directory=os.path.relpath(inputs, output),
        batch_plan=os.path.relpath(args.plan.resolve(), output)), indent=2) + '\n')
    start = time.monotonic()
    status = dict(state='running', start_utc=datetime.now(timezone.utc).isoformat(), seeds={})

    def save():
        status['elapsed_s'] = time.monotonic() - start
        tmp = output / 'batch_status.tmp'
        tmp.write_text(json.dumps(status, indent=2) + '\n')
        tmp.replace(output / 'batch_status.json')

    save()
    try:
        with (output / 'batch.log').open('w', buffering=1) as log:
            lock = threading.Lock()

            def run_seed(seed):
                run = output / f'seed_{seed:05d}'
                completed = subprocess.run([sys.executable, '-B', str(inputs / 'run.py'), str(run),
                    '--inputs', str(inputs), '--seed', str(seed)], stdout=log, stderr=subprocess.STDOUT)
                with lock:
                    status['seeds'][str(seed)] = dict(returncode=completed.returncode,
                        state=json.loads((run / 'status.json').read_text())['state'] if (run / 'status.json').exists() else 'failed')
                    save()

            with ThreadPoolExecutor(args.workers or plan.get('workers', 1)) as pool:
                list(pool.map(run_seed, plan['seeds']))
            status['state'] = 'summarizing'
            save()
            extra = ['--seed1', str(repo / plan['seed1'])] if plan.get('seed1') else []
            subprocess.run([sys.executable, '-B', str(inputs / 'summarize_batch.py'), str(output), *extra],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
            status['state'] = 'complete' if all(r['state'] == 'complete' for r in status['seeds'].values()) else 'completed_with_failures'
    except BaseException as error:
        status.update(state='failed', error=repr(error))
        raise
    finally:
        status['end_utc'] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
