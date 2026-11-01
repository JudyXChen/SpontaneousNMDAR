"""Compare original and calibrated NMDAR block kinetics after a test pulse.

This is an exact population calculation for the ten-state receptor, not a
spatial MCell simulation. It isolates the effect of changing k_B(V, M) and
k_U(V) while keeping every other transition rate unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.nmdar.compare_conductance import (  # noqa: E402
    block_rate as calibrated_block_rate,
    load_rates,
    unblock_rate as calibrated_unblock_rate,
)
from scripts.plot_style import apply_style, condition_color, figure, save, save_legend  # noqa: E402

STATES = ('C0', 'C1', 'C2', 'D', 'O', 'B', 'C2B', 'DB', 'C1B', 'C0B')
INDEX = {state: index for index, state in enumerate(STATES)}
ELEMENTARY_CHARGE_C = 1.602176634e-19
MODELED_CA_CONDUCTANCE_S = 4.5e-15
MODELED_REVERSAL_MV = 3.0


def original_block_rate(voltage_mV: float, magnesium_mM: float) -> float:
    return 1200.0 * magnesium_mM * math.exp(-voltage_mV / 17.0)


def original_unblock_rate(voltage_mV: float) -> float:
    return 10800.0 * math.exp(voltage_mV / 47.0)


def generator(voltage_mV: float, magnesium_mM: float, glutamate_M: float,
              rates: dict[str, float], calibrated: bool) -> np.ndarray:
    q = np.zeros((len(STATES), len(STATES)), dtype=float)

    def edge(source: str, target: str, rate: float) -> None:
        source_index = INDEX[source]
        target_index = INDEX[target]
        q[target_index, source_index] += rate
        q[source_index, source_index] -= rate

    def pair(a: str, b: str, forward: float, reverse: float) -> None:
        edge(a, b, forward)
        edge(b, a, reverse)

    a1 = rates['k_C0C1'] * glutamate_M
    a2 = rates['k_C1C2'] * glutamate_M
    pair('C0', 'C1', a1, rates['k_C1C0'])
    pair('C1', 'C2', a2, rates['k_C2C1'])
    pair('C2', 'D', rates['k_C2D'], rates['k_DC2'])
    pair('C2', 'O', rates['k_C20'], rates['k_0C2'])
    if calibrated:
        block = calibrated_block_rate(voltage_mV, magnesium_mM)
        unblock = calibrated_unblock_rate(voltage_mV)
    else:
        block = original_block_rate(voltage_mV, magnesium_mM)
        unblock = original_unblock_rate(voltage_mV)
    pair('O', 'B', block, unblock)
    pair('C2B', 'B', rates['k_C20b'], rates['k_0C2b'])
    pair('C2B', 'DB', rates['k_C2D'], rates['k_DC2'])
    pair('C1B', 'C2B', a2, rates['k_C2C1'])
    pair('C0B', 'C1B', a1, rates['k_C1C0'])
    return q


def propagate(q: np.ndarray, initial: np.ndarray,
              elapsed_s: np.ndarray) -> np.ndarray:
    """Evaluate exp(q*t) @ initial for all nonnegative elapsed times."""
    eigenvalues, eigenvectors = np.linalg.eig(q)
    inverse = np.linalg.inv(eigenvectors)
    coefficients = inverse @ initial
    result = np.empty((len(elapsed_s), len(initial)), dtype=float)
    for index, elapsed in enumerate(elapsed_s):
        state = eigenvectors @ (np.exp(eigenvalues * elapsed) * coefficients)
        state = np.real_if_close(state, tol=1000).real
        state[np.abs(state) < 1e-14] = 0.0
        result[index] = state / state.sum()
    if np.min(result) < -1e-10:
        raise ArithmeticError('Negative state probability from propagation')
    return result


def pulse_response(voltage_mV: float, magnesium_mM: float,
                   glutamate_uM: float, pulse_ms: float, duration_ms: float,
                   step_ms: float, calibrated: bool,
                   rates: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    time_ms = np.arange(0.0, duration_ms + step_ms / 2.0, step_ms)
    pulse_end = int(round(pulse_ms / step_ms))
    if not np.isclose(time_ms[pulse_end], pulse_ms):
        raise ValueError('Pulse duration must be an integer multiple of the output step')
    initial = np.zeros(len(STATES))
    initial[INDEX['C0']] = 1.0
    q_on = generator(voltage_mV, magnesium_mM, glutamate_uM * 1e-6,
                     rates, calibrated)
    during = propagate(q_on, initial, time_ms[:pulse_end + 1] * 1e-3)
    q_off = generator(voltage_mV, magnesium_mM, 0.0, rates, calibrated)
    after_elapsed = (time_ms[pulse_end + 1:] - pulse_ms) * 1e-3
    after = propagate(q_off, during[-1], after_elapsed)
    return time_ms, np.vstack((during, after))


def summarize(time_ms: np.ndarray, probabilities: np.ndarray,
              receptor_count: int, voltage_mV: float) -> dict[str, float]:
    open_probability = probabilities[:, INDEX['O']]
    expected_open = receptor_count * open_probability
    probability_any_open = 1.0 - np.power(1.0 - open_probability, receptor_count)
    integrated_open_receptor_s = float(
        np.trapz(expected_open, time_ms * 1e-3))
    event_rate_per_open = (MODELED_CA_CONDUCTANCE_S
                           * (MODELED_REVERSAL_MV - voltage_mV) * 1e-3
                           / (2.0 * ELEMENTARY_CHARGE_C))
    return {
        'peak_single_receptor_open_probability': float(open_probability.max()),
        'peak_expected_open_receptors': float(expected_open.max()),
        'peak_probability_any_receptor_open': float(probability_any_open.max()),
        'integrated_open_receptor_ms': integrated_open_receptor_s * 1e3,
        'expected_modeled_ca_events': integrated_open_receptor_s * event_rate_per_open,
    }


def rate_diagnostics(voltage_mV: float, magnesium_mM: float,
                     rates: dict[str, float], calibrated: bool) -> dict[str, float]:
    if calibrated:
        block = calibrated_block_rate(voltage_mV, magnesium_mM)
        unblock = calibrated_unblock_rate(voltage_mV)
    else:
        block = original_block_rate(voltage_mV, magnesium_mM)
        unblock = original_unblock_rate(voltage_mV)
    open_exit = block + rates['k_0C2']
    blocked_exit = unblock + rates['k_0C2b']
    return {
        'k_B_per_s': block,
        'k_U_per_s': unblock,
        'mean_open_state_sojourn_us': 1e6 / open_exit,
        'mean_blocked_B_state_sojourn_us': 1e6 / blocked_exit,
        'probability_B_exits_by_unblocking': unblock / blocked_exit,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--voltage-mv', type=float, default=-67.0)
    parser.add_argument('--magnesium-mm', type=float, default=1.0)
    parser.add_argument('--glutamate-um', type=float, default=100.0)
    parser.add_argument('--pulse-ms', type=float, default=1.0)
    parser.add_argument('--duration-ms', type=float, default=35.0)
    parser.add_argument('--step-ms', type=float, default=0.01)
    parser.add_argument('--receptors', type=int, default=10)
    parser.add_argument(
        '--output-dir', type=Path,
        default=(ROOT / 'archive/NMDAR_old/comparison/'
                 'kinetics/100uM_Glu_1ms'))
    args = parser.parse_args()
    if not (args.magnesium_mm >= 0 and args.glutamate_um > 0
            and 0 < args.pulse_ms <= args.duration_ms and args.step_ms > 0
            and args.receptors > 0):
        raise ValueError('Invalid test protocol')

    rates = load_rates()
    responses = {}
    summaries = {}
    kinetic_rates = {}
    for name, calibrated in (('original', False), ('calibrated', True)):
        time_ms, probabilities = pulse_response(
            args.voltage_mv, args.magnesium_mm, args.glutamate_um,
            args.pulse_ms, args.duration_ms, args.step_ms, calibrated, rates)
        responses[name] = probabilities
        summaries[name] = summarize(
            time_ms, probabilities, args.receptors, args.voltage_mv)
        kinetic_rates[name] = rate_diagnostics(
            args.voltage_mv, args.magnesium_mm, rates, calibrated)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    protocol = {
        'voltage_mV': args.voltage_mv,
        'magnesium_mM': args.magnesium_mm,
        'glutamate_uM': args.glutamate_um,
        'pulse_ms': args.pulse_ms,
        'duration_ms': args.duration_ms,
        'step_ms': args.step_ms,
        'receptors': args.receptors,
        'initial_state': 'C0',
    }
    payload = {
        'protocol': protocol,
        'kinetic_rates': kinetic_rates,
        'summaries': summaries,
    }
    (args.output_dir / 'summary.json').write_text(
        json.dumps(payload, indent=2) + '\n')

    with (args.output_dir / 'open_probability.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow([
            'time_ms', 'p_open_original', 'p_open_calibrated',
            'expected_open_10_original', 'expected_open_10_calibrated'])
        for index, time in enumerate(time_ms):
            original = responses['original'][index, INDEX['O']]
            calibrated = responses['calibrated'][index, INDEX['O']]
            writer.writerow([
                f'{time:.8g}', f'{original:.12g}', f'{calibrated:.12g}',
                f'{args.receptors * original:.12g}',
                f'{args.receptors * calibrated:.12g}'])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    apply_style()
    trace_colors = {
        'original': condition_color('blue', 'bright'),
        'calibrated': condition_color('blue', 'dark'),
    }
    fig, ax = figure()
    for name in ('original', 'calibrated'):
        expected_open = args.receptors * responses[name][:, INDEX['O']]
        ax.plot(time_ms, expected_open, color=trace_colors[name], label=name.capitalize())
    ax.set(xlabel='Time (ms)', ylabel='Expected open NMDAR')
    save(fig, args.output_dir / 'expected_open.png', legend_name=None)

    handles = [Line2D([0], [0], color=trace_colors[name], label=name.capitalize())
               for name in ('original', 'calibrated')]
    save_legend(handles, args.output_dir / 'legend.png', figsize=(3, 0.65), ncol=2)
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
