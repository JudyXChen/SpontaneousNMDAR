"""Resolve immutable inputs separately from simulation outputs."""
import hashlib
import json
from pathlib import Path


def inputs_for(run):
    run = Path(run).resolve()
    reference = run / 'input_reference.json'
    if reference.exists():
        return (run / json.loads(reference.read_text())['input_directory']).resolve()
    return run  # Historical layout, before the approved migration.


def verify_inputs(inputs):
    for name, expected in json.loads((inputs / 'input_hashes.json').read_text()).items():
        actual = hashlib.sha256((inputs / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed snapshot input: ' + name)


def repository():
    return next(p for p in Path(__file__).resolve().parents if (p / 'scripts/plot_style.py').exists())
