"""Find where a model condition's recorded output lives under results/."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def results_for(model):
    relative = Path(model).resolve().relative_to(ROOT / 'model')
    return ROOT / 'results' / relative.parent / 'conditions' / relative.name
