"""Resolve active and archived models independently of recorded outputs."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / 'archive/NMDAR_old'


def model_family(name):
    return ARCHIVE / 'model' if name == 'NMDAR_old' else ROOT / 'model' / name


def result_family(name):
    return ARCHIVE / 'results' if name == 'NMDAR_old' else ROOT / 'results' / name


def model_for(key):
    family, relative = str(key).split('/', 1)
    return model_family(family) / relative


def results_for(model):
    path = Path(model).resolve()
    for source, output in [(ROOT / 'model', ROOT / 'results'),
                           (ARCHIVE / 'model', ARCHIVE / 'results')]:
        if path.is_relative_to(source):
            relative = path.relative_to(source)
            return output / relative.parent / 'conditions' / relative.name
    return path
