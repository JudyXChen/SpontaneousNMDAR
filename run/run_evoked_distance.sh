#!/usr/bin/env bash
# Experiments 2 and 3: calibrated evoked distance comparisons, three seeds each.
# Requires MCell4 and its compatible Python (set MCELL_PATH and PYTHON).
# Usage: bash run_evoked_distance.sh [--smoke] [--jobs 3] [--output NEW_DIRECTORY]
# Full option list: bash run_evoked_distance.sh --help
# Outputs are isolated under results/calcium_total/evoked/distance_seeds2_4/<unique run>/ by default.
# This file embeds a standard-library Python scheduler; no extra packages needed.
set -euo pipefail
TASK_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
# Use the verified local CellBlender runtime when no explicit runtime is set.
# Server installations can still supply MCELL_PATH and PYTHON in the environment.
LOCAL_RESOURCES=/Applications/Blender-4.4.3-CellBlender/Blender.app/Contents/Resources/4.4
if [[ -z "${MCELL_PATH:-}" && -f "$LOCAL_RESOURCES/extensions/system/cellblender/extensions/mcell/lib/mcell.so" ]]; then
    export MCELL_PATH="$LOCAL_RESOURCES/extensions/system/cellblender/extensions/mcell"
fi
if [[ -z "${PYTHON:-}" && "${MCELL_PATH:-}" == "$LOCAL_RESOURCES/extensions/system/cellblender/extensions/mcell" ]]; then
    export PYTHON="$LOCAL_RESOURCES/python/bin/python3.11"
fi
exec "${PYTHON:-python3}" -B -u - "$TASK_ROOT" "$@" <<'PY'
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import tempfile

ROOT = Path(sys.argv[1])
parser = argparse.ArgumentParser(
    prog='run_evoked_distance.sh',
    description="Run calibrated evoked distance experiments with at most five concurrent MCell processes.",
    epilog="Set MCELL_PATH to the directory containing lib/mcell.so and PYTHON to its compatible interpreter. "
           "Each invocation creates fresh outputs; existing output directories are refused. "
           "Model defaults are preserved; timestep is fixed at 1e-7 s and glutamate sites at both.")
parser.add_argument('--experiment', choices=('2', '3', 'both'), default='both',
                    help='2: 10 NMDAR/cluster; 3: 3 NMDAR/cluster (default: both)')
parser.add_argument('--seeds', nargs='+', type=int, default=[2, 3, 4],
                    help='random seeds (default: 2 3 4)')
parser.add_argument('--jobs', type=int, default=3, help='concurrent simulations, 1 to 5 (default: 3)')
parser.add_argument('--iterations', type=int, help='multiple of 100 (default: 350000 = 35 ms)')
parser.add_argument('--smoke', action='store_true', help='1000 iterations = 0.1 ms for every selected case/seed')
parser.add_argument('--input-root', type=Path, help='snapshot parent; default: model/shared/snapshots')
parser.add_argument('--output', type=Path, help='new output directory; default is timestamped under results/')
args = parser.parse_args(sys.argv[2:])
if not 1 <= args.jobs <= 5:
    parser.error('--jobs must be between 1 and 5')
if not args.seeds or len(set(args.seeds)) != len(args.seeds) or any(s < 1 or s > 2147483647 for s in args.seeds):
    parser.error('seeds must be unique positive integers, at most 2147483647')
if args.smoke and args.iterations is not None:
    parser.error('use either --smoke or --iterations')
iterations = args.iterations if args.iterations is not None else (1000 if args.smoke else 350000)
if iterations < 100 or iterations % 100:
    parser.error('--iterations must be a positive multiple of 100 (observable sampling interval)')
dt = 1e-7
cases = []
if args.experiment in ('2', 'both'):
    cases += [f'2_cluster_{d}nm' for d in (50, 100, 200)]
if args.experiment in ('3', 'both'):
    cases += [f'2_cluster_{d}nm_3_NMDAR_each' for d in (50, 100, 200)]

raw_mcell_path = os.environ.get('MCELL_PATH')
if not raw_mcell_path:
    parser.error('set MCELL_PATH to your MCell4 installation directory')
mcell_path = Path(raw_mcell_path).expanduser().resolve()
if not (mcell_path / 'lib/mcell.so').is_file():
    parser.error(f'missing {mcell_path}/lib/mcell.so')
env = os.environ.copy()
# Ignore inherited experiment overrides: the command-line options are authoritative.
for key in list(env):
    if key.startswith('MCELL_'):
        del env[key]
env.update(MCELL_PATH=str(mcell_path), MCELL_ITERATIONS=str(iterations),
           MCELL_TIME_STEP=str(dt), MCELL_GLU_SITES='both', PYTHONDONTWRITEBYTECODE='1',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
           NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
probe = subprocess.run([
    sys.executable, '-B', '-c',
    "import os,sys,json; sys.path.insert(0,os.path.join(os.environ['MCELL_PATH'],'lib')); "
    "import mcell; print(json.dumps({'version':str(mcell.__version__),'module':mcell.__file__}))"
], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if probe.returncode:
    parser.error('MCell import failed with this PYTHON; check its Python ABI and shared libraries:\n' + probe.stderr)

source_files = {}
for case in cases:
    source = ROOT / 'model/calcium_total/evoked' / case
    required = ['model.py', 'model.bngl', 'parameters.py', 'customization.py',
                'instantiation.py', 'observables.py', 'subsystem.py', 'shared.py', 'geometry/__init__.py']
    for name in required:
        if not (source / name).is_file():
            parser.error(f'missing model input: {source / name}')
    source_files[case] = sorted(list(source.glob('*.py')) + list(source.glob('*.bngl'))
                                + list((source / 'geometry').glob('*.py')))

stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
default_output = (Path(tempfile.mkdtemp(prefix='nmdar-smoke-')) / 'run'
                  if args.smoke else ROOT / 'results/calcium_total/evoked/distance_seeds2_4' / ('run_' + stamp))
out = (args.output or default_output).resolve()
if out.exists():
    parser.error(f'output directory already exists; choose a NEW directory: {out}')
out.mkdir(parents=True)
default_inputs = (Path(tempfile.mkdtemp(prefix='nmdar-inputs-'))
                  if args.smoke else ROOT / 'model/shared/snapshots')
inputs = ((args.input_root or default_inputs) / stamp).resolve()
if inputs.is_relative_to(ROOT / 'results'):
    parser.error('Input snapshots must be outside results/')
inputs.mkdir(parents=True, exist_ok=False)

def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)

def git_info(*arguments):
    result = subprocess.run(['git', '-C', str(ROOT), *arguments], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return result.stdout.strip() if result.returncode == 0 else 'unavailable'

# Snapshot only model inputs, never old outputs, bytecode, or checkpoints.
hashes = {}
for case, paths in source_files.items():
    source = ROOT / 'model/calcium_total/evoked' / case
    for path in paths:
        relative = path.relative_to(source)
        target = inputs / case / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / 'model/shared/legacy/model.py' if relative.as_posix() == 'model.py' else path, target)
        hashes[f'{case}/{relative}'] = hashlib.sha256(target.read_bytes()).hexdigest()
shutil.copy2(ROOT / 'run/run_evoked_distance.sh', inputs / 'run_evoked_distance.sh')
write_json(out / 'input_reference.json', {'input_directory': os.path.relpath(inputs, out)})
manifest = dict(created_utc=stamp, git_commit=git_info('rev-parse', 'HEAD'),
                git_status=git_info('status', '--short'), python=sys.version,
                python_executable=sys.executable, mcell=probe.stdout.strip(),
                cases=cases, seeds=args.seeds, jobs=args.jobs, iterations=iterations,
                time_step_s=dt, duration_s=iterations * dt, glutamate_sites='both',
                smoke=args.smoke, input_sha256=hashes,
                command=sys.argv[2:], output=str(out))
write_json(out / 'manifest.json', manifest)
records = [dict(case=case, seed=seed, status='pending') for case in cases for seed in args.seeds]

def save_status():
    write_json(out / 'status.json', records)
    lines = ['case\tseed\tstatus\texit_code\telapsed_s']
    for record in records:
        lines.append('\t'.join(str(record.get(k, '')) for k in
                               ('case', 'seed', 'status', 'exit_code', 'elapsed_s')))
    (out / 'status.tsv').write_text('\n'.join(lines) + '\n')

def validate_counts(work, seed):
    # The three observables used by both existing distance-comparison scripts.
    # Verify the entire sampling grid, not just the process exit code.
    for name in ('NMDAR_World.dat', 'NMDAR_Ca_flux_World.dat', 'Ca2_Mushroom_transparent.dat'):
        path = work / 'react_data' / f'seed_{seed:05d}' / name
        rows = path.read_text().splitlines()
        if len(rows) != iterations // 100 + 1:
            raise ValueError(f'{name}: expected {iterations // 100 + 1} samples, got {len(rows)}')
        for i, line in enumerate(rows):
            t, value = map(float, line.split())
            if not math.isfinite(value) or not math.isclose(t, i * 100 * dt, rel_tol=0, abs_tol=1e-10):
                raise ValueError(f'{name}: invalid sample at row {i + 1}')

active = {}
pending = iter(records)
exhausted = False
interrupted = False

def stop(signum, frame):
    raise KeyboardInterrupt

signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)
save_status()
print(f'Output: {out}\nCases: {len(cases)}; seeds: {args.seeds}; total runs: {len(records)}; workers: {args.jobs}')
print(f'Duration: {iterations * dt * 1000:g} ms; timestep: {dt:g} s; glutamate sites: both')
print('Each simulation has its own directory. Progress: status.tsv; details: CASE/seed_NNNNN/run.log')
try:
    while active or not exhausted:
        while len(active) < args.jobs and not exhausted:
            record = next(pending, None)
            if record is None:
                exhausted = True
                break
            work = out / record['case'] / f"seed_{record['seed']:05d}"
            work.mkdir(parents=True)
            write_json(work / 'input_reference.json', {'input_directory': os.path.relpath(inputs / record['case'], work)})
            log = (work / 'run.log').open('w')
            run_env = dict(env, MCELL_SEED=str(record['seed']))
            try:
                process = subprocess.Popen([sys.executable, '-B', '-u', str(inputs / record['case'] / 'model.py'), '-seed', str(record['seed'])],
                                           cwd=work, env=run_env, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True)
            finally:
                log.close()
            record.update(status='running', pid=process.pid, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
            active[process] = (record, work, time.monotonic())
            print(f"START {record['case']} seed {record['seed']} (PID {process.pid})")
            save_status()
        for process, (record, work, started) in list(active.items()):
            code = process.poll()
            if code is None:
                continue
            record.update(exit_code=code, elapsed_s=round(time.monotonic() - started, 2))
            record['status'] = 'failed'
            if code == 0:
                try:
                    validate_counts(work, record['seed'])
                    record['status'] = 'complete'
                except (OSError, ValueError) as error:
                    record['error'] = str(error)
            else:
                record['error'] = f'MCell exited with status {code}; see run.log'
            print(f"{record['status'].upper()} {record['case']} seed {record['seed']} ({record['elapsed_s']} s)")
            del active[process]
            save_status()
        if active:
            time.sleep(0.2)
except KeyboardInterrupt:
    interrupted = True
    print('\nStopping active simulations; partial outputs will be retained.')
finally:
    # Also clean up children if an unexpected scheduler error occurs.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    for process in active:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    for process, (record, work, started) in active.items():
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
        record.update(status='interrupted', exit_code=process.returncode,
                      elapsed_s=round(time.monotonic() - started, 2))
    for record in records:
        if record['status'] == 'pending':
            record['status'] = 'not_run'
    save_status()
completed = sum(r['status'] == 'complete' for r in records)
print(f'Completed {completed}/{len(records)} runs. Status: {out / "status.tsv"}')
sys.exit(130 if interrupted else (0 if completed == len(records) else 1))
PY
