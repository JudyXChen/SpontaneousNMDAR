"""Find where a model condition's recorded output lives under results/."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def results_for(model):
    """results/<family>/<protocol>/runs/<case>: holds seed_NNNNN/ run folders and plots/."""
    relative = Path(model).resolve().relative_to(ROOT / 'model')
    return ROOT / 'results' / relative.parent / 'runs' / relative.name
