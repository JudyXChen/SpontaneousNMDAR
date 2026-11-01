"""Build and run experiments 2, 3-4, or 5-6 with a bounded process queue."""
import argparse
import datetime
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import tempfile

from build import COMMON, MODEL_ROOT, load_condition, materialize
from validate_run import validate_run

ROOT = MODEL_ROOT.parent


def suite_leaves(suite):
    calibrated = MODEL_ROOT / 'calcium_total'
    if suite == 'spont_1_vs_2':
        return [calibrated / 'spont/1_cluster_10_NMDAR', calibrated / 'spont/2_cluster_100nm']
    family = 'spont'
    return [calibrated / family / (f'2_cluster_{d}nm' + ('_3_NMDAR_each' if n == 3 else ''))
            for n in (10, 3) for d in (50, 100, 200)]


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def main(condition=None):
    parser = argparse.ArgumentParser(description=__doc__)
    if condition is None:
        parser.add_argument('--suite', required=True, choices=('spont_1_vs_2', 'spont_distance'))
    parser.add_argument('--seeds', nargs='+', type=int, help='default: seed 1 for experiment 2; seeds 1-4 for sweeps')
    parser.add_argument('--jobs', type=int, default=3, help='concurrent simulations (default: 3)')
    parser.add_argument('--smoke', action='store_true', help='compress timeline to 0.1 ms baseline + 0.1 ms response')
    parser.add_argument('--prepare-only', action='store_true', help='write self-contained inputs without running MCell')
    parser.add_argument('--list', action='store_true', help='print protocol and conditions without writing files')
    parser.add_argument('--input-root', type=Path, help='snapshot parent; default: model/shared/snapshots')
    parser.add_argument('--output', type=Path, help='new, nonexisting output directory')
    args = parser.parse_args()
    suite = args.suite if condition is None else 'single_condition'
    seeds = args.seeds if args.seeds is not None else ([1] if suite in ('spont_1_vs_2', 'single_condition') else [1, 2, 3, 4])
    if args.jobs < 1 or args.jobs > 64:
        parser.error('--jobs must be between 1 and 64')
    if len(set(seeds)) != len(seeds) or any(s < 1 or s > 2147483647 for s in seeds):
        parser.error('seeds must be unique positive integers <= 2147483647')
    leaves = [condition] if condition is not None else suite_leaves(suite)
    configs = [load_condition(leaf, smoke=args.smoke) for leaf in leaves]
    print(f'Suite: {suite}; {len(configs)} conditions; seeds {seeds}; workers {args.jobs}', flush=True)
    for config in configs:
        print(f"  {config['condition']}: {config['duration_s'] * 1000:g} ms; "
              f"NMDAR {[c['population'] + ':' + str(c['receptors']) for c in config['clusters']]}; "
              f"VGCC {config['vgcc_populations']}", flush=True)
    if args.list:
        return
    env = {k: v for k, v in os.environ.items() if not k.startswith('MCELL_')}
    env.update(PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    runtime = {'python': sys.version, 'executable': sys.executable}
    if not args.prepare_only:
        raw = os.environ.get('MCELL_PATH', '')
        if not raw:
            parser.error('MCELL_PATH is not set; use a supplied shell launcher or configure the MCell runtime')
        mcell_path = Path(raw).expanduser().resolve()
        if not (mcell_path / 'lib/mcell.so').is_file():
            parser.error(f'missing {mcell_path}/lib/mcell.so')
        env['MCELL_PATH'] = str(mcell_path)
        probe = subprocess.run([sys.executable, '-B', '-c',
                                "import sys,os;sys.path.insert(0,os.path.join(os.environ['MCELL_PATH'],'lib'));"
                                "import mcell;print(mcell.__version__);print(mcell.__file__)"],
                               env=env, capture_output=True, text=True)
        if probe.returncode:
            parser.error('MCell import failed; use its compatible Python:\n' + probe.stderr)
        runtime['mcell'] = probe.stdout.strip()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    default_output = (Path(tempfile.mkdtemp(prefix='nmdar-smoke-')) / 'run'
                      if args.smoke else ROOT / 'results' / suite / ('run_' + stamp))
    out = (args.output or default_output).resolve()
    if out.exists():
        parser.error('Output directory already exists: ' + str(out))
    out.mkdir(parents=True)
    default_inputs = (Path(tempfile.mkdtemp(prefix='nmdar-inputs-'))
                      if args.smoke else COMMON / 'snapshots')
    inputs = ((args.input_root or default_inputs) / stamp).resolve()
    if inputs.is_relative_to(ROOT / 'results'):
        parser.error('Input snapshots must be outside results/')
    inputs.mkdir(parents=True, exist_ok=False)
    write_json(out / 'input_reference.json', {'input_directory': os.path.relpath(inputs, out)})
    records, hashes = [], {}
    for config in configs:
        case = Path(config['condition']).name
        hashes[case] = materialize(config, inputs / case)
        records.extend(dict(case=case, seed=seed, status='pending') for seed in seeds)
    revision = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    write_json(out / 'manifest.json', dict(suite=suite, created_utc=stamp, seeds=seeds, jobs=args.jobs,
                                         smoke=args.smoke, configs=configs, runtime=runtime,
                                         git_revision=revision.stdout.strip(), input_sha256=hashes))
    # Capture the orchestration and validation code as well as all simulation inputs.
    (inputs / 'launcher').mkdir()
    for path in COMMON.glob('*.py'):
        shutil.copy2(path, inputs / 'launcher' / path.name)
    print('Output: ' + str(out), flush=True)

    def save_status():
        write_json(out / 'status.json', records)
        fields = ('case', 'seed', 'status', 'exit_code', 'elapsed_s')
        (out / 'status.tsv').write_text('\t'.join(fields) + '\n' + ''.join(
            '\t'.join(str(row.get(k, '')) for k in fields) + '\n' for row in records))

    if args.prepare_only:
        for row in records:
            row['status'] = 'prepared'
        save_status()
        return
    config_by_case = {Path(c['condition']).name: c for c in configs}
    active, pending, exhausted, interrupted = {}, iter(records), False, False

    def stop(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    save_status()
    try:
        while active or not exhausted:
            while len(active) < args.jobs and not exhausted:
                row = next(pending, None)
                if row is None:
                    exhausted = True
                    break
                work = out / row['case'] / f"seed_{row['seed']:05d}"
                work.mkdir(parents=True)
                write_json(work / 'input_reference.json', {'input_directory': os.path.relpath(inputs / row['case'], work)})
                with (work / 'run.log').open('w') as log:
                    process = subprocess.Popen([sys.executable, '-B', '-u', str(inputs / row['case'] / 'model.py'), '-seed', str(row['seed'])],
                                               cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
                                               start_new_session=True)
                row.update(status='running', pid=process.pid)
                active[process] = (row, work, time.monotonic())
                print(f"START {row['case']} seed {row['seed']}", flush=True)
                save_status()
            for process, (row, work, started) in list(active.items()):
                code = process.poll()
                if code is None:
                    continue
                row.update(status='failed', exit_code=code, elapsed_s=round(time.monotonic() - started, 2))
                if code == 0:
                    try:
                        validate_run(work, row['seed'], config_by_case[row['case']])
                        row['status'] = 'complete'
                    except (OSError, ValueError) as error:
                        row['error'] = str(error)
                else:
                    row['error'] = 'MCell failed; inspect run.log'
                print(f"{row['status'].upper()} {row['case']} seed {row['seed']} ({row['elapsed_s']} s)", flush=True)
                del active[process]
                save_status()
            if active:
                time.sleep(.2)
    except KeyboardInterrupt:
        interrupted = True
        print('Stopping active simulations; partial outputs are retained.', flush=True)
    finally:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        for process in active:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        for process, (row, work, started) in active.items():
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            row.update(status='interrupted', exit_code=process.returncode)
        for row in records:
            if row['status'] == 'pending':
                row['status'] = 'not_run'
        save_status()
    completed = sum(row['status'] == 'complete' for row in records)
    print(f'Completed {completed}/{len(records)} runs. Status: {out / "status.tsv"}', flush=True)
    raise SystemExit(130 if interrupted else (0 if completed == len(records) else 1))


if __name__ == '__main__':
    main()
