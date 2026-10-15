#!/usr/bin/env python3
"""Build and run this condition using the shared model code."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'shared'))
from batch import main
if __name__ == '__main__':
    main(condition=Path(__file__).resolve().parent)
