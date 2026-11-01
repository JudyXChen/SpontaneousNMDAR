"""Matched 35 ms NMDAR cluster conditions and their output paths."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'archive/NMDAR_old/results/evoked'
ALL_COMPARISON = 'compare_3cases'
PAIR_COMPARISON = 'compare_10NMDAR'
CASES = (
    ('NMDAR_old/evoked/1_cluster_10_NMDAR',
     'One cluster, 10 NMDAR', (1681, 1686, 1688)),
    ('NMDAR_old/evoked/1_cluster_20_NMDAR',
     'One cluster, 20 NMDAR',
     (1681, 1686, 1688, 1690, 1696, 1719)),
    ('NMDAR_old/evoked/2_cluster_100nm',
     'Two clusters, 10 + 10 NMDAR',
     ((1681, 1686, 1688), (1687, 1691, 1694))),
)
