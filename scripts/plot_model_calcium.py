"""Plot calcium traces from the IP3R homog MCell run, averaged across seeds.

Produces figures in the matching model branch under results/:
  ca_{er,cyto}.png -- one compartment per file
  open_{ryr,ip3r}.png -- one channel per file
  ip3r_gate.png -- IP3R open vs single-subunit-in-110
  er_ode.png -- MCell vs ODE (with --ode-csv)
                                    Only written when --ode-csv is given.

Data files (per seed, columns = [time_s, molecule_count]):
  Ca2_ER_SA.dat        -> ER luminal calcium
  Ca2_World.dat        -> cytosolic free calcium
  RyR_O_World.dat      -> open RyR channels
  IP3R_open_World.dat  -> open IP3R channels (all three subunits in 110)
  IP3R_s1_110_World.dat-> channels with subunit s1 in 110 (old-model analogue)

Usage:
  python scripts/plot_model_calcium.py MODEL_DIR [--ode-csv PATH]

The ODE overlay is skipped unless --ode-csv points at an ERPM traces.csv; it lives
in a different repo and is not present on a compute server.
"""
import argparse
import csv
import glob
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for
from scripts.plot_style import condition_output_dir, apply_style, colors, figure, legend_below, panel_path, save, shade
apply_style()

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument('model_dir', type=Path, help='Model leaf; recorded data is resolved under results/')
ap.add_argument('--output-dir', type=Path, default=None)
ap.add_argument(
    "--ode-csv", default=None, metavar="PATH",
    help="ERPM ODE traces.csv for the MCell-vs-ODE ER overlay. Lives in a "
         "different repo; the overlay is skipped when this is not given.")
ARGS = ap.parse_args()
MODEL_DIR = ARGS.model_dir.resolve()
REACT = str(results_for(MODEL_DIR) / 'react_data')
if ARGS.output_dir is None:
    OUTDIR = str(condition_output_dir(MODEL_DIR))
else:
    OUTDIR = str(ARGS.output_dir.resolve())
os.makedirs(OUTDIR, exist_ok=True)


def load_mean(fname):
    """Return (time_s, mean, std, n_seeds) for `fname` averaged across seeds."""
    traces, time = [], None
    for sd in sorted(glob.glob(os.path.join(REACT, "seed_*"))):
        p = os.path.join(sd, fname)
        if not os.path.exists(p):
            continue
        d = np.loadtxt(p)
        if d.ndim != 2 or d.shape[1] != 2 or len(d) < 2:
            continue  # A seed with only its initial state is incomplete.
        if time is None:
            time = d[:, 0]
        elif not np.array_equal(time, d[:, 0]):
            raise ValueError(f"Seed time grid differs: {p}")
        traces.append(d[:, 1])
    if not traces:
        raise FileNotFoundError(f"No complete traces for {fname} in {REACT}")
    traces = np.array(traces)
    return time, traces.mean(axis=0), traces.std(axis=0), len(traces)


def load_ode(csv_path, column):
    """Return (time_s, values) for one column of the ODE traces CSV."""
    with open(csv_path) as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    ci, ti = header.index(column), header.index("t")
    t = np.array([float(r[ti]) for r in rows[1:]])
    v = np.array([float(r[ci]) for r in rows[1:]])
    return t, v


def _split_traces(specs, out):
    """Save each mean +/-1 SD trace in its own figure."""
    for suffix, (tms, mean, std, color, ylabel) in specs:
        fig, ax = figure()
        ax.plot(tms, mean, color=color, label="_nolegend_")
        ax.fill_between(tms, mean - std, mean + std, color=color, alpha=0.2,
                        label="_nolegend_")
        ax.set_xlabel("time (ms)")
        ax.set_ylabel(ylabel)
        save(fig, panel_path(os.path.join(OUTDIR, out), suffix))


def plot_er_cyto():
    """ER and cytosolic calcium (molecule counts), mean +/-1 SD across seeds."""
    t, er_mean, er_std, n = load_mean("Ca2_ER_SA.dat")
    _, cy_mean, cy_std, _ = load_mean("Ca2_World.dat")
    tms = t * 1e3  # s -> ms

    print(f"ER   Ca2 ({n} seeds): {er_mean[0]:.0f} -> {er_mean[-1]:.0f} molecules")
    print(f"Cyto Ca2 ({n} seeds): {cy_mean[0]:.0f} -> {cy_mean[-1]:.0f} molecules")

    _split_traces(
        [("er", (tms, er_mean, er_std, colors(2)[0], "ER Ca$^{2+}$ (molecules)")),
         ("cyto", (tms, cy_mean, cy_std, colors(2)[1], "Cytosolic Ca$^{2+}$ (molecules)"))],
        "ca.png")


def plot_open_channels():
    """RyR and IP3R open counts vs time, mean +/-1 SD across seeds."""
    t, ryr_mean, ryr_std, n = load_mean("RyR_O_World.dat")
    _, ip3r_mean, ip3r_std, _ = load_mean("IP3R_open_World.dat")
    tms = t * 1e3  # s -> ms

    print(f"RyR  open ({n} seeds): peak {ryr_mean.max():.1f} channels")
    print(f"IP3R open ({n} seeds): peak {ip3r_mean.max():.1f} channels")

    _split_traces(
        [("ryr", (tms, ryr_mean, ryr_std, colors(2)[0], "Open RyR (count)")),
         ("ip3r", (tms, ip3r_mean, ip3r_std, colors(2)[1], "Open IP3R (count)"))],
        "open.png")


def plot_ip3r_open():
    """Open IP3R channels vs single-subunit-in-110, showing the effect of cubing.

    IP3R(s1~110) counts channels with subunit s1 in 110 regardless of the other two
    -- the direct analogue of the published model's IP3R_110 observable, which let a
    single subunit conduct. IP3R(s1~110,s2~110,s3~110) requires all three, which is
    the DYK gate. The gap between the two curves is the bug being fixed.
    """
    t, open_mean, open_std, n = load_mean("IP3R_open_World.dat")
    _, s1_mean, s1_std, _ = load_mean("IP3R_s1_110_World.dat")
    tms = t * 1e3

    subs = {}
    for s in ("s1", "s2", "s3"):
        try:
            _, m_, _, _ = load_mean(f"IP3R_{s}_110_World.dat")
            subs[s] = m_
        except (IndexError, ValueError):
            pass

    print(f"IP3R open channels ({n} seeds): peak {open_mean.max():.3f}, "
          f"final {open_mean[-1]:.3f}")
    print(f"IP3R single subunit in 110:     peak {s1_mean.max():.3f}, "
          f"final {s1_mean[-1]:.3f}")
    if len(subs) == 3:
        # Statistical-equivalence check: the three subunits are independent and
        # identically parameterized, so these agree in distribution but not exactly.
        means = {s: v.mean() for s, v in subs.items()}
        spread = max(means.values()) - min(means.values())
        print(f"subunit time-means s1/s2/s3: "
              + " ".join(f"{s}={v:.3f}" for s, v in means.items())
              + f"  (spread {spread:.3f})")

    fig, ax = figure()
    one_color, open_color = shade("blue"), shade("blue", "dark")
    ax.plot(tms, s1_mean, color=one_color,
            label="1 subunit in 110\n(published gate)")
    ax.fill_between(tms, s1_mean - s1_std, s1_mean + s1_std, color=one_color, alpha=0.2)
    ax.plot(tms, open_mean, color=open_color,
            label="3 subunits in 110 (DYK gate)")
    ax.fill_between(tms, open_mean - open_std, open_mean + open_std,
                    color=open_color, alpha=0.2)
    ax.axhline(1.0, color="k", lw=0.6, ls=":", alpha=0.6)
    ax.text(tms[-1], 1.0, " 1 channel", va="bottom", ha="right",
            alpha=0.7)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("channels (count)")
    legend_below(ax)
    ax.margins(x=0)
    path = os.path.join(OUTDIR, "ip3r_gate.png")
    save(fig, path)


def plot_er_ode_vs_mcell(ode_csv):
    """MCell ER depletion vs ODE (SOCE off), each normalized to its t=0 value."""
    if not ode_csv:
        print("Skipping ODE overlay: no --ode-csv given")
        return
    if not os.path.exists(ode_csv):
        print(f"Skipping ODE overlay: {ode_csv} not found")
        return
    ODE_CSV = ode_csv

    t, er_mean, er_std, _ = load_mean("Ca2_ER_SA.dat")
    mc_tms = t * 1e3
    mc_norm = er_mean / er_mean[0]
    mc_norm_std = er_std / er_mean[0]

    ode_t, ode_caer = load_ode(ODE_CSV, "Ca_ER")
    ode_tms = ode_t * 1e3
    ode_norm = ode_caer / ode_caer[0]

    print(f"MCell mean ER: {er_mean[0]:.0f} -> {er_mean[-1]:.0f} "
          f"({100 * (1 - mc_norm[-1]):.1f}% depletion)")
    print(f"ODE  Ca_ER:    {ode_caer[0] * 1e6:.1f} uM -> {ode_caer[-1] * 1e6:.1f} uM "
          f"({100 * (1 - ode_norm[-1]):.1f}% depletion)")

    fig, ax = figure()
    mcell_color, ode_color = shade("blue", "dark"), shade("blue")
    ax.plot(mc_tms, mc_norm, color=mcell_color, label="MCell mean")
    ax.fill_between(mc_tms, mc_norm - mc_norm_std, mc_norm + mc_norm_std,
                    color=mcell_color, alpha=0.2, label="MCell ±1 SD")
    ax.plot(ode_tms, ode_norm, color=ode_color, ls="--", label="ODE (SOCE off)")
    ax.axhline(1.0, color="k", lw=0.5, alpha=0.4)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("Normalized ER Ca$^{2+}$")
    legend_below(ax)

    path = os.path.join(OUTDIR, "er_ode.png")
    save(fig, path)


if __name__ == "__main__":
    plot_er_cyto()
    plot_open_channels()
    plot_ip3r_open()
    plot_er_ode_vs_mcell(ARGS.ode_csv)
