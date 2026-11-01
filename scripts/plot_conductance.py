"""Compare calibrated ten-state NMDAR conductance with Jahr-Stevens Eq. 5.

This is a stationary, fixed-glutamate calculation, not an MCell trajectory or
a new conductance equation. Only O conducts; B remains nonconducting and can
reversibly unblock to O. The two block/unblock rates were calibrated at 100 uM
glutamate to match Jahr-Stevens Eq. 5 for 0.1 and 1.25 mM Mg. Because the
stationary curve identifies only k_B/k_U, the fit preserves the original
geometric mean sqrt(k_B*k_U) at every voltage and divides the log-rate change
equally between the two rates. The solid curve is P_O(V, Mg) / P_O(V, 0), not
raw P_O. This comparison does not assume a calcium-current fraction.

References:
  Vargas-Caballero & Robinson (2004), doi:10.1523/JNEUROSCI.1380-04.2004
  Bartol et al. (2015), doi:10.3389/fnsyn.2015.00017
  Jahr & Stevens (1990), doi:10.1523/JNEUROSCI.10-09-03178.1990
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_style import apply_style, figure, shade, save
MODEL = ROOT / 'model/calcium_total/evoked/1_cluster_10_NMDAR'
STATES = ('C0', 'C1', 'C2', 'D', 'O', 'B', 'C2B', 'DB', 'C1B', 'C0B')
STATE_INDEX = {name: i for i, name in enumerate(STATES)}
PARAMETERS = ('k_C0C1', 'k_C1C0', 'k_C1C2', 'k_C2C1', 'k_C20',
              'k_0C2', 'k_C20b', 'k_0C2b', 'k_C2D', 'k_DC2')

# Calibrated at 100 uM glutamate against Jahr-Stevens Eq. 5. The stationary
# ten-state model has
#   P_O(M) / P_O(0) = 1 / (1 + K * k_B(V, M) / k_U(V)),
# where K=2.8398079368402547 is the blocked/unblocked branch-weight ratio at
# 100 uM glutamate. Matching Eq. 5 fixes k_B/k_U but not the absolute rates.
# These values impose the additional minimal-change constraint
# k_B(new)*k_U(new) = k_B(original)*k_U(original) at every voltage.
BLOCK_PREFACTOR_PER_MM = 1130.6388608604038
BLOCK_VOLTAGE_SCALE_MV = 20.091025673263097
UNBLOCK_PREFACTOR = 11462.546042453892
UNBLOCK_VOLTAGE_SCALE_MV = 81.78933360630568


def load_rates(path: Path = MODEL / 'model.bngl') -> dict[str, float]:
    rates = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) == 2 and words[0] in PARAMETERS:
            rates[words[0]] = float(words[1])
    missing = set(PARAMETERS) - rates.keys()
    if missing:
        raise ValueError(f'Missing original NMDAR rates: {sorted(missing)}')
    return rates


def validate_calibrated_voltage_tables(path: Path = MODEL / 'parameters.py') -> None:
    """Ensure the MCell tables implement the calibrated 1 mM rate pair."""
    found = {}
    for node in ast.parse(path.read_text()).body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in ('k_B_rate', 'k_U_rate')):
            found[node.targets[0].id] = np.asarray(ast.literal_eval(node.value), dtype=float)
    if set(found) != {'k_B_rate', 'k_U_rate'}:
        raise ValueError('Calibrated block/unblock tables were not found')
    block, unblock = found['k_B_rate'], found['k_U_rate']
    if block.shape != unblock.shape or block.ndim != 2 or block.shape[1] != 2:
        raise ValueError('Calibrated block/unblock tables have incompatible shapes')
    if not np.array_equal(block[:, 0], unblock[:, 0]):
        raise ValueError('Calibrated block/unblock tables use different times')
    voltage = -BLOCK_VOLTAGE_SCALE_MV * np.log(
        block[:, 1] / BLOCK_PREFACTOR_PER_MM)
    expected_unblock = UNBLOCK_PREFACTOR * np.exp(
        voltage / UNBLOCK_VOLTAGE_SCALE_MV)
    if not np.allclose(unblock[:, 1], expected_unblock, rtol=1e-12, atol=1e-9):
        raise ValueError('Tables do not match the calibrated voltage rates')


def block_rate(voltage_mV: float, magnesium_mM: float) -> float:
    return (BLOCK_PREFACTOR_PER_MM * magnesium_mM
            * math.exp(-voltage_mV / BLOCK_VOLTAGE_SCALE_MV))


def unblock_rate(voltage_mV: float) -> float:
    return UNBLOCK_PREFACTOR * math.exp(
        voltage_mV / UNBLOCK_VOLTAGE_SCALE_MV)


def transition_matrix(voltage_mV: float, magnesium_mM: float,
                      glutamate_M: float, rates: dict[str, float]) -> np.ndarray:
    """Column-source generator for the full ten-state stochastic scheme."""
    if not (math.isfinite(voltage_mV) and math.isfinite(magnesium_mM)
            and math.isfinite(glutamate_M)):
        raise ValueError('Voltage, magnesium, and glutamate must be finite')
    if magnesium_mM < 0 or glutamate_M <= 0:
        raise ValueError('Magnesium must be nonnegative and glutamate positive')
    q = np.zeros((len(STATES), len(STATES)), dtype=float)

    def edge(source: str, target: str, rate: float) -> None:
        i, j = STATE_INDEX[source], STATE_INDEX[target]
        q[j, i] += rate
        q[i, i] -= rate

    def pair(a: str, b: str, forward: float, reverse: float) -> None:
        edge(a, b, forward)
        edge(b, a, reverse)

    a1 = rates['k_C0C1'] * glutamate_M
    a2 = rates['k_C1C2'] * glutamate_M
    d1, d2 = rates['k_C1C0'], rates['k_C2C1']
    pair('C0', 'C1', a1, d1)
    pair('C1', 'C2', a2, d2)
    pair('C2', 'D', rates['k_C2D'], rates['k_DC2'])
    pair('C2', 'O', rates['k_C20'], rates['k_0C2'])
    pair('O', 'B', block_rate(voltage_mV, magnesium_mM),
         unblock_rate(voltage_mV))
    pair('C2B', 'B', rates['k_C20b'], rates['k_0C2b'])
    pair('C2B', 'DB', rates['k_C2D'], rates['k_DC2'])
    pair('C1B', 'C2B', a2, d2)
    pair('C0B', 'C1B', a1, d1)
    return q


def stationary_open_probability(voltage_mV: float, magnesium_mM: float,
                                glutamate_M: float, rates: dict[str, float]) -> float:
    q = transition_matrix(voltage_mV, magnesium_mM, glutamate_M, rates)
    balance = q.copy()
    balance[-1, :] = 1.0
    target = np.zeros(len(STATES))
    target[-1] = 1.0
    p = np.linalg.solve(balance, target)
    if np.any(p < -1e-10) or not np.isclose(p.sum(), 1.0, atol=1e-10):
        raise ArithmeticError('Invalid stationary occupancy')
    return float(p[STATE_INDEX['O']])


def jahr_factor(voltage_mV: np.ndarray, magnesium_mM: float) -> np.ndarray:
    """Jahr-Stevens high-Mg approximation, with V in mV and Mg in mM."""
    return 1.0 / (1.0 + (magnesium_mM / 3.57) * np.exp(-0.062 * voltage_mV))


def make_comparison(glutamate_uM: float, magnesium_mM: list[float],
                    voltage_min: float, voltage_max: float, points: int,
                    output_dir: Path) -> None:
    if not math.isfinite(glutamate_uM) or glutamate_uM <= 0:
        raise ValueError('Glutamate must be positive and finite')
    if (not magnesium_mM or len(magnesium_mM) > 3
            or len(set(magnesium_mM)) != len(magnesium_mM)
            or any(not math.isfinite(m) or m <= 0 for m in magnesium_mM)):
        raise ValueError('Use one to three distinct, positive, finite Mg concentrations')
    if not (math.isfinite(voltage_min) and math.isfinite(voltage_max)
            and voltage_min < voltage_max and points >= 3):
        raise ValueError('Invalid voltage grid')
    validate_calibrated_voltage_tables()
    rates = load_rates()
    glutamate_M = glutamate_uM * 1e-6
    voltage = np.linspace(voltage_min, voltage_max, points)
    # With Mg=0, only block rates depend on V, so this baseline is V-independent.
    open_no_magnesium = stationary_open_probability(0.0, 0.0, glutamate_M, rates)
    curves = {}
    for magnesium in magnesium_mM:
        label = f'{magnesium:g} mM'
        kinetic = np.array([stationary_open_probability(v, magnesium, glutamate_M, rates)
                            / open_no_magnesium for v in voltage])
        jahr = jahr_factor(voltage, magnesium)
        curves[label] = {'ten_state': kinetic, 'jahr': jahr}

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / 'comparison.csv'
    column_names = {'ten_state': 'po_ratio', 'jahr': 'jahr_g'}
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['voltage_mV'] + [f'{column_names[model]}_{m:g}mM'
                                          for m in magnesium_mM for model in ('ten_state', 'jahr')])
        for i, v in enumerate(voltage):
            writer.writerow([f'{v:.8g}'] + [f'{curves[f"{m:g} mM"][model][i]:.10g}'
                                                for m in magnesium_mM
                                                for model in ('ten_state', 'jahr')])

    plot_comparison(voltage, curves, output_dir)


def plot_saved_comparison(output_dir: Path) -> None:
    """Replot recorded conductance values without changing either kinetic model."""
    data = np.genfromtxt(output_dir / 'comparison.csv', delimiter=',', names=True, deletechars='')
    curves = {}
    for column in data.dtype.names:
        if column.startswith('po_ratio_'):
            magnesium = column.removeprefix('po_ratio_').removesuffix('mM')
            curves[f'{float(magnesium):g} mM'] = {
                'ten_state': data[column], 'jahr': data[f'jahr_g_{magnesium}mM']}
    plot_comparison(data['voltage_mV'], curves, output_dir)


def plot_comparison(voltage, curves, output_dir: Path) -> None:
    """Render either saved original or newly calculated calibrated curves."""
    magnesium_mM = [float(label.split()[0]) for label in curves]
    voltage_min, voltage_max = float(voltage[0]), float(voltage[-1])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    apply_style()
    fig, ax = figure()
    # One hue for Mg; higher concentration gets a darker shade.
    tones = ('main', 'dark', 'deep')
    mg_colors = {magnesium: shade('blue', tone) for magnesium, tone in
                 zip(sorted(magnesium_mM), tones[-len(magnesium_mM):])}
    for magnesium in sorted(magnesium_mM):
        color = mg_colors[magnesium]
        row = curves[f'{magnesium:g} mM']
        ax.plot(voltage, row['ten_state'], color=color)
        ax.plot(voltage, row['jahr'], color=color,
                linestyle=(0, (5, 3)))
    ax.set(xlim=(voltage_min, voltage_max), ylim=(-0.02, 1.03),
           xlabel='Membrane voltage (mV)', ylabel='Normalized conductance')
    ax.set_yticks(np.linspace(0, 1, 6))
    style_handles = [
        Line2D([0], [0], color='black',
               linestyle=(0, (5, 3)),
               label=r'Jahr-Stevens 1990: $g(V)$'),
        Line2D([0], [0], color='black',
               label=r'MCell: $P_O/P_O(M=0)$')]
    mg_handles = [Line2D([0], [0], color=color,
                         label=f'{magnesium:g} mM Mg' + r'$^{2+}$')
                  for magnesium, color in sorted(mg_colors.items())]
    # Place the two legend groups below the full-size panel.
    fig.legend(handles=style_handles, loc='upper left',
               bbox_to_anchor=(0.02, -0.02), frameon=False,
               handlelength=2, handletextpad=0.4, borderaxespad=0)
    fig.legend(handles=mg_handles, loc='upper left',
               bbox_to_anchor=(0.72, -0.02), frameon=False,
               handlelength=1, handletextpad=0.4, borderaxespad=0)
    save(fig, output_dir / 'comparison.png', legend_name='legend.png')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--from-csv', action='store_true',
                        help='Replot comparison.csv without recalculating kinetics.')
    parser.add_argument('--glutamate-um', type=float, default=100.0)
    parser.add_argument('--magnesium-mm', type=float, nargs='+', default=[1.25, 10.0, 0.1])
    parser.add_argument('--voltage-min', type=float, default=-100.0)
    parser.add_argument('--voltage-max', type=float, default=60.0)
    parser.add_argument('--points', type=int, default=401)
    parser.add_argument(
        '--output-dir', type=Path,
        default=ROOT / 'results/calcium_total/conductance')
    args = parser.parse_args()
    if args.from_csv:
        plot_saved_comparison(args.output_dir)
    else:
        make_comparison(args.glutamate_um, args.magnesium_mm,
                        args.voltage_min, args.voltage_max,
                        args.points, args.output_dir)
    print('Wrote title-free PNG figures and CSV to', args.output_dir)
