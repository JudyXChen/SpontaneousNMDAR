"""Run an exported evoked condition using separate input and output directories."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

COMMON = Path(__file__).resolve().parent
ROOT = COMMON.parent.parent


def materialize(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    files = sorted(list(source.glob('*.py')) + list(source.glob('*.bngl'))
                   + list((source / 'geometry').glob('*.py')))
    if not (source / 'model.bngl').exists():
        raise ValueError('Expected an exported evoked model leaf with model.bngl')
    for path in files:
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        runtime = source.parents[2] / 'shared/legacy/model.py'
        shutil.copy2(runtime if relative.as_posix() == 'model.py' else path, target)
    return {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in destination.rglob('*') if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path)
    parser.add_argument('output', type=Path, help='new result directory; existing directories are refused')
    parser.add_argument('seeds', type=int, nargs='*', default=[1])
    args = parser.parse_args()
    seeds = args.seeds or [1]
    if any(s < 1 for s in seeds) or len(set(seeds)) != len(seeds):
        parser.error('Seeds must be unique positive integers')
    out = args.output.resolve()
    if out.exists():
        parser.error('Output already exists; choose a new directory')
    if out.is_relative_to(ROOT / 'model'):
        parser.error('Use results/ for output, not model/')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    inputs = COMMON / 'snapshots' / stamp / args.model.name
    hashes = materialize(args.model, inputs)
    out.mkdir(parents=True)
    (out / 'input_reference.json').write_text(json.dumps(dict(
        input_directory=os.path.relpath(inputs, out), source=str(args.model),
        input_sha256=hashes, seeds=seeds), indent=2) + '\n')
    (out / 'logs').mkdir()
    for seed in seeds:
        with (out / 'logs' / f'seed_{seed}.log').open('w') as log:
            subprocess.run([sys.executable, '-B', str(inputs / 'model.py'), '-seed', str(seed)],
                           cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(out)


if __name__ == '__main__':
    main()
