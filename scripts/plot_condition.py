"""Plot one condition's recorded data (calcium_total).

  plot_condition.py traces MODEL_DIR   PSD layout and calcium/NMDAR traces
  plot_condition.py calcium MODEL_DIR  ER and cytosolic calcium, open RyR and IP3R
  plot_condition.py flux MODEL_DIR     calcium fluxes and compartment budgets

Figures go to results/ for the model, or --output-dir. Calcium and flux average
all seeds. Flux values are rebuilt from count observables and the model rate laws,
in ions per second, positive into the cytosol.
"""
from __future__ import annotations

import argparse
import ast
import csv
import glob
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for
from scripts.plot_style import (condition_output_dir, apply_style, colors, condition_color,
                                figure, legend_below, panel_path, save, shade)


# ---- traces ----

PANELS = (
    ('NMDAR_Ca_flux_World.dat',
     'Cumulative modeled NMDAR\nCa$^{2+}$ entry (firings)', 'ca_entry'),
    ('NMDAR_World.dat', 'Open NMDAR (count)', 'open_nmdar'),
    ('Ca2_Mushroom_transparent.dat',
     'Spine cytosolic free Ca$^{2+}$ (ions)', 'spine_cytosolic_ca'),
    ('Ca2_head_neck.dat', 'Head and neck free Ca$^{2+}$ (ions)', 'ca_head_neck'),
    ('Ca2_Spine_base.dat', 'Spine-base free Ca$^{2+}$ (ions)', 'ca_base'),
)

TRACE_COLORS = {
    'calcium_total/evoked/1_cluster_10_NMDAR': condition_color('blue', 'dark'),
    'calcium_total/evoked/2_cluster_50nm': condition_color('purple', 'dark'),
    'calcium_total/evoked/2_cluster_100nm': condition_color('orange', 'dark'),
    'calcium_total/evoked/2_cluster_200nm': condition_color('blue', 'dark'),
    'calcium_total/evoked/2_cluster_50nm_3_NMDAR_each':
        condition_color('purple', 'dark'),
    'calcium_total/evoked/2_cluster_100nm_3_NMDAR_each':
        condition_color('orange', 'dark'),
    'calcium_total/evoked/2_cluster_200nm_3_NMDAR_each':
        condition_color('blue', 'dark'),
}


def trace_color(model_dir: Path) -> str:
    """Return the approved condition color, or the default for other cases."""
    try:
        key = model_dir.relative_to(ROOT / 'model').as_posix()
    except ValueError:
        return colors(1)[0]
    return TRACE_COLORS.get(key, colors(1)[0])


def literal_assignments(path: Path, names: set[str]) -> dict[str, object]:
    values: dict[str, object] = {}
    for node in ast.parse(path.read_text()).body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in names):
            values[node.targets[0].id] = ast.literal_eval(node.value)
    missing = names - values.keys()
    if missing:
        raise ValueError(f'Missing assignments in {path}: {sorted(missing)}')
    return values


def load_trace(model_dir: Path, seed: int, file_name: str) -> np.ndarray:
    path = results_for(model_dir) / 'react_data' / f'seed_{seed:05d}' / file_name
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f'Expected a two-column trace: {path}')
    if len(values) != 3501 or not np.isclose(values[-1, 0], 0.035, atol=1e-10):
        raise ValueError(f'Expected a complete 35 ms trace with 3501 samples: {path}')
    return values


def positions(model_dir: Path, seed: int, iteration: int, species: str) -> np.ndarray:
    path = (results_for(model_dir) / 'viz_data' / f'seed_{seed:05d}'
            / f'Scene.ascii.{iteration:06d}.dat')
    if iteration == 0 and not path.exists():
        path = path.with_name('Scene.ascii.0000.dat')
    points = []
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0] == species:
                points.append([float(value) for value in fields[2:5]])
    return np.asarray(points, dtype=float).reshape(-1, 3)


def plot_layout(model_dir: Path, output_dir: Path, seed: int) -> None:
    geometry = literal_assignments(
        model_dir / 'geometry/Mushroom_dendrite.py',
        {'Mushroom_dendrite_vertex_list', 'Mushroom_dendrite_wall_list',
         'Mushroom_dendrite_PSD_wall_indices'},
    )
    layout = literal_assignments(
        model_dir / 'instantiation.py', {'CLUSTER_WALL_INDICES'})
    vertices = np.asarray(geometry['Mushroom_dendrite_vertex_list'])
    walls = geometry['Mushroom_dendrite_wall_list']
    regions = layout['CLUSTER_WALL_INDICES']
    if regions and isinstance(regions[0], int):
        regions = (regions,)

    def polygons(indices) -> np.ndarray:
        return np.asarray([vertices[walls[index], :2] * 1000 for index in indices])

    centers = np.asarray([polygons(indices).mean(axis=(0, 1)) for indices in regions])
    receptors = positions(model_dir, seed, 0, 'NMDAR_0')[:, :2] * 1000
    assigned = np.argmin(
        np.linalg.norm(receptors[:, None, :] - centers[None, :, :], axis=2), axis=1)
    palette = colors(3)
    fig, ax = figure()
    ax.add_collection(PolyCollection(
        polygons(geometry['Mushroom_dendrite_PSD_wall_indices']),
        facecolor='#eeeeee', edgecolor='#d2d2d2', linewidth=0.3))
    for number, (indices, center, color) in enumerate(
            zip(regions, centers, palette[:2]), 1):
        selected = assigned == number - 1
        ax.add_collection(PolyCollection(
            polygons(indices), facecolor=color, alpha=0.2,
            edgecolor=color, linewidth=1))
        ax.scatter(receptors[selected, 0], receptors[selected, 1], s=24,
                   color=color, label=f'Cluster {number} (n={selected.sum()})',
                   zorder=3)
        ax.scatter(center[0], center[1], s=90, marker='x', color=palette[2],
                   linewidths=2,
                   label='Glutamate site' if number == 1 else None, zorder=4)
    if len(centers) == 2:
        separation = np.linalg.norm(centers[0] - centers[1])
        print(f'Realized cluster-center separation: {separation:.3f} nm')
    ax.set_aspect('equal')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    legend_below(ax)
    save(fig, output_dir / 'psd_layout.png')


def run_traces(args):
    model_dir = args.model_dir.resolve()
    output_dir = condition_output_dir(model_dir) if args.output_dir is None else args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_layout(model_dir, output_dir, args.seed)
    color = trace_color(model_dir)
    for file_name, ylabel, suffix in PANELS:
        values = load_trace(model_dir, args.seed, file_name)
        fig, ax = figure()
        ax.plot(values[:, 0] * 1000, values[:, 1], color=color, drawstyle='steps-post')
        ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(0, 35))
        ax.set_xticks((0, 10, 20, 30))
        if file_name == 'NMDAR_World.dat':
            ax.set_ylim(bottom=-0.05)
        save(fig, output_dir / f'{suffix}.png')
        print(f'{file_name}: final={values[-1, 1]:g}, peak={values[:, 1].max():g}')


# ---- calcium ----

def run_calcium(ARGS):
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


    plot_er_cyto()
    plot_open_channels()
    plot_ip3r_open()
    plot_er_ode_vs_mcell(ARGS.ode_csv)


# ---- flux ----

def run_flux(ARGS):
    N_A = 6.02214076e23          # molecules / mol

    Ca_ER_IC = 150e-6            # M   (Release_Site_12)
    Ca_cyto_IC = 1e-7            # M   (Release_Site_1)

    k_PMCA_3 = 12.0              # 1/s   PMCA1 -> PMCA           (extrusion)
    k_PMCA_leak = 4.3           # 1/s   PMCA  -> PMCA + Ca       (leak)
    k_NCX_3 = 600.0             # 1/s   NCX1  -> NCX             (extrusion)
    k_NCX_leak = 19.4           # 1/s   NCX   -> NCX + Ca        (leak)
    k_SERCA_leak = 0.1608       # 1/s   Ca_ER + SA_mem -> Ca     (ER leak)
    k_FluxRyR = 1.09e9          # M^-1 s^-1  RyR_O Ca flux
    k_IP3R_Ca_flux = 1.19e8     # M^-1 s^-1  IP3R_110 Ca flux
    k_Ca_N = 50.0               # 1/s   Ca2 -> 0                 (cytosolic sink)

    # Total pumps = placement density (per um^2) x membrane area (um^2).
    # PMCA & NCX sit on head+neck+cylinder_receptors of Mushroom_dendrite; the summed
    # triangle area of those regions is 5.617 um^2 (computed from the geometry mesh).
    PM_AREA = 5.617                  # um^2
    N_PMCA_TOTAL = 998.0 * PM_AREA   # ~5606  (PMCA density 998 /um^2)
    N_NCX_TOTAL = 142.0 * PM_AREA    # ~798   (NCX  density 142 /um^2)

    MODEL_DIR = ARGS.model_dir.resolve()
    REACT = str(results_for(MODEL_DIR) / 'react_data')
    PARAMS = str(MODEL_DIR / 'parameters.py')
    if ARGS.output_dir is None:
        OUTDIR = str(condition_output_dir(MODEL_DIR))
    else:
        OUTDIR = str(ARGS.output_dir.resolve())
    os.makedirs(OUTDIR, exist_ok=True)


    # --------------------------------------------------------------------------
    # loaders / helpers
    # --------------------------------------------------------------------------
    def load_seed_matrix(fname):
        """Return (time_s, M) where M has shape (n_seeds, n_time) for `fname`."""
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
        return time, np.array(traces)


    def parse_rate_table(name):
        """Extract a `name = [[t, rate], ...]` variable-rate table from parameters.py."""
        src = open(PARAMS).read()
        m = re.search(re.escape(name) + r"\s*=\s*\[", src)
        i = m.end() - 1
        depth = 0
        for j in range(i, len(src)):
            if src[j] == "[":
                depth += 1
            elif src[j] == "]":
                depth -= 1
                if depth == 0:
                    return np.array(eval(src[i:j + 1]), dtype=float)
        raise RuntimeError("could not parse rate table " + name)


    def smooth(x, win=201):
        """Odd-window moving average with edge padding (for the noisy ER derivative)."""
        if len(x) <= 3:
            return x
        win = min(win, (len(x) // 2) * 2 - 1)
        if win < 3:
            return x
        pad = win // 2
        xp = np.pad(x, pad, mode="edge")
        return np.convolve(xp, np.ones(win) / win, mode="valid")


    def cumtrapz(y, x):
        """Cumulative trapezoidal integral of y over x, starting at 0."""
        out = np.zeros_like(y, dtype=float)
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
        return out


    # --------------------------------------------------------------------------
    # load counts (per seed) and derive volumes
    # --------------------------------------------------------------------------
    t, ca_er = load_seed_matrix("Ca2_ER_SA.dat")
    _, ca_cyto = load_seed_matrix("Ca2_World.dat")
    _, vscc = load_seed_matrix("VSCC_World.dat")
    _, nmdar = load_seed_matrix("NMDAR_World.dat")
    _, pmca1 = load_seed_matrix("PMCA1_World.dat")
    _, ncx1 = load_seed_matrix("NCX1_World.dat")
    _, ryr_o = load_seed_matrix("RyR_O_World.dat")
    _, ip3r_open = load_seed_matrix("IP3R_open_World.dat")
    _, serca_count = load_seed_matrix("SERCA_flux_World.dat")
    _, serca_leak_count = load_seed_matrix("SERCA_leak_World.dat")
    _, ryr_count = load_seed_matrix("RyR_flux_net_World.dat")
    _, ip3r_count = load_seed_matrix("IP3R_flux_World.dat")
    _, cabf = load_seed_matrix("CaBf_World.dat")
    _, cabm = load_seed_matrix("CaBm_Mushroom_transparent.dat")

    n_seeds = ca_er.shape[0]
    tms = t * 1e3  # s -> ms

    # effective compartment volumes implied by the initial conditions (litres)
    V_ER = ca_er[:, 0].mean() / (N_A * Ca_ER_IC)
    V_CYTO = ca_cyto[:, 0].mean() / (N_A * Ca_cyto_IC)
    conc_er = ca_er / (N_A * V_ER)          # [M], per seed
    conc_cyto = ca_cyto / (N_A * V_CYTO)    # [M], per seed

    # time-dependent single-channel rates (np.interp holds the last tabulated value
    # past the 0-50 ms stimulus window, matching MCell's variable-rate behaviour)
    vscc_tab = parse_rate_table("vscc_rate_scaled")
    nmdar_tab = parse_rate_table("nmdar_rate_scaled")
    vscc_rate = np.interp(t, vscc_tab[:, 0], vscc_tab[:, 1])
    nmdar_rate = np.interp(t, nmdar_tab[:, 0], nmdar_tab[:, 1])

    # --------------------------------------------------------------------------
    # per-seed fluxes in the cytosol-positive convention (+ into cyto, - out of cyto)
    # --------------------------------------------------------------------------
    def count_rate(counts):
        """Convert cumulative reaction counts to per-seed firing rates (ions/s)."""
        return np.vstack([np.gradient(row, t) for row in counts])


    # Directly measured reaction-count fluxes. Positive magnitudes point ER -> cytosol,
    # except SERCA uptake, whose positive magnitude points cytosol -> ER.
    mag_serca = count_rate(serca_count)
    mag_serca_l = count_rate(serca_leak_count)
    mag_ryr = count_rate(ryr_count)
    mag_ip3r = count_rate(ip3r_count)

    # Independent comparisons: the old IP3R rate-law reconstruction and ER
    # mass-balance SERCA residual. Neither is used as the measured flux.
    mag_ip3r_reconstructed = k_IP3R_Ca_flux * ip3r_open * conc_er
    dNdt_er = np.vstack([np.gradient(ca_er[s], t) for s in range(n_seeds)])
    mag_serca_residual = dNdt_er + mag_serca_l + mag_ryr + mag_ip3r

    J = {
        "VSCC":    vscc_rate[None, :] * vscc,                    # ECS -> cyto  (+)
        "NMDAR":   nmdar_rate[None, :] * nmdar,                  # ECS -> cyto  (+)
        "PMCA_L":  k_PMCA_leak * (N_PMCA_TOTAL - pmca1),         # ECS -> cyto  (+)
        "NCX_L":   k_NCX_leak * (N_NCX_TOTAL - ncx1),            # ECS -> cyto  (+)
        "PMCA":   -k_PMCA_3 * pmca1,                             # cyto -> ECS  (-)
        "NCX":    -k_NCX_3 * ncx1,                               # cyto -> ECS  (-)
        "SERCA_L": mag_serca_l,                                  # ER  -> cyto  (+)
        "RyR":     mag_ryr,                                      # ER  -> cyto  (+)
        "IP3R":    mag_ip3r,                                     # ER  -> cyto  (+)
        "SERCA":  -mag_serca,                                    # cyto -> ER   (-)
        "SERCA_RESIDUAL": -mag_serca_residual,                   # cross-check only
    }

    # --------------------------------------------------------------------------
    # figure 1: one signed panel per flux
    # --------------------------------------------------------------------------
    IN = "ECS → cytosol"      # plasma-membrane influx  (+)
    OUT = "cytosol → ECS"     # plasma-membrane efflux  (-)
    REL = "ER → cytosol"      # ER release              (+)
    UP = "cytosol → ER"       # ER uptake               (-)
    CAT_COLOR = {IN: shade("blue"), OUT: shade("purple"),
                 REL: shade("orange"), UP: shade("blue", "dark")}

    panels = [
        ("$J_{VSCC}$",    "VSCC",   IN),
        ("$J_{NMDAR}$",   "NMDAR",  IN),
        ("$J_{PMCA,L}$",  "PMCA_L", IN),
        ("$J_{NCX,L}$",   "NCX_L",  IN),
        ("$J_{PMCA}$",    "PMCA",   OUT),
        ("$J_{NCX}$",     "NCX",    OUT),
        ("$J_{SERCA,L}$", "SERCA_L", REL),
        ("$J_{RyR}$",     "RyR",    REL),
        ("$J_{IP3R}$",    "IP3R",   REL),
        ("$J_{SERCA}$",   "SERCA",  UP),
        ("$J_{SERCA}$ residual", "SERCA_RESIDUAL", UP),
    ]


    def plot_fluxes():
        print(f"Fluxes from {n_seeds} seeds "
              f"(V_ER={V_ER*1e15:.4f} um^3, V_cyto={V_CYTO*1e15:.4f} um^3):")
        for name, key, cat in panels:
            fig, ax = figure()
            color = CAT_COLOR[cat]
            m, s = J[key].mean(axis=0), J[key].std(axis=0)
            print(f"  {name:14s} {cat:16s} signed peak = {m[np.argmax(np.abs(m))]:+.3e} ions/s")
            ax.plot(tms, m, color=color, lw=1.6)
            ax.fill_between(tms, m - s, m + s, color=color, alpha=0.2, lw=0)
            ax.axhline(0, color="k", lw=0.6, alpha=0.4)
            ax.margins(x=0)
            ax.set_ylabel(f"{name} flux (ions s$^{{-1}}$)")
            ax.set_xlabel("time (ms)")
            save(fig, panel_path(os.path.join(OUTDIR, "flux.png"), key.lower()))


    def plot_ip3r_comparison():
        """Direct IP3R reaction-count flux versus the rate-law reconstruction."""
        direct = mag_ip3r.mean(axis=0)
        reconstructed = mag_ip3r_reconstructed.mean(axis=0)
        fig, ax = figure()
        ax.plot(tms, direct, color=shade("blue", "dark"), label="Direct count")
        ax.plot(tms, reconstructed, color=shade("blue"), ls="--",
                label="Rate-law reconstruction")
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("IP3R Ca$^{2+}$ flux (ions s$^{-1}$)")
        legend_below(ax)
        ax.margins(x=0)
        out = os.path.join(OUTDIR, "ip3r_flux.png")
        save(fig, out)


    # --------------------------------------------------------------------------
    # figure 2: cumulative Ca budget per compartment (mean across seeds)
    # --------------------------------------------------------------------------
    def plot_budget():
        """Decompose the measured Ca change in each compartment into cumulative fluxes."""
        Jm = {k: v.mean(axis=0) for k, v in J.items()}   # mean signed flux
        cy = ca_cyto.mean(axis=0)
        er = ca_er.mean(axis=0)

        figE, axE = figure()
        figC, axC = figure()

        # ---- ER lumen: exact balance of SERCA (in) vs leak/RyR/IP3R (out) ----
        # ER gains -J_SERCA (cyto->ER) and loses the release fluxes (ER->cyto).
        er_terms = [
            ("SERCA uptake", -cumtrapz(Jm["SERCA"], t),   shade("blue", "dark")),
            ("SERCA leak",   -cumtrapz(Jm["SERCA_L"], t), shade("blue")),
            ("RyR",          -cumtrapz(Jm["RyR"], t),     shade("orange")),
            ("IP3R",         -cumtrapz(Jm["IP3R"], t),    shade("purple")),
        ]
        er_sum = sum(v for _, v, _ in er_terms)
        er_meas = er - er[0]
        for label, v, c in er_terms:
            axE.plot(tms, v, color=c, lw=1.6, label=label)
        axE.plot(tms, er_meas, color="k", label="Measured")
        axE.plot(tms, er_sum, color="0.5", ls="--", label="Flux sum")
        axE.axhline(0, color="k", lw=0.6, alpha=0.4)
        axE.set_xlabel("time (ms)")
        axE.set_ylabel("ER cumulative Ca (ions)")
        legend_below(axE, ncol=2)
        axE.margins(x=0)

        # ---- cytosol: membrane fluxes + buffering + decay, with a residual ----
        buffer_net = cabf.mean(axis=0) + cabm.mean(axis=0)      # Ca bound to buffers
        cy_terms = [
            ("PM influx",
             cumtrapz(Jm["VSCC"] + Jm["NMDAR"] + Jm["PMCA_L"] + Jm["NCX_L"], t), shade("blue")),
            ("PM efflux",
             cumtrapz(Jm["PMCA"] + Jm["NCX"], t), shade("purple")),
            ("ER release",
             cumtrapz(Jm["SERCA_L"] + Jm["RyR"] + Jm["IP3R"], t), shade("orange")),
            ("ER uptake",
             cumtrapz(Jm["SERCA"], t), shade("blue", "dark")),
            ("Buffering",
             -(buffer_net - buffer_net[0]), shade("orange", "dark")),
            ("Ca decay",
             -cumtrapz(k_Ca_N * cy, t), shade("purple", "dark")),
        ]
        cy_sum = sum(v for _, v, _ in cy_terms)
        cy_meas = cy - cy[0]
        residual = cy_meas - cy_sum
        for label, v, c in cy_terms:
            axC.plot(tms, v, color=c, lw=1.6, label=label)
        axC.plot(tms, cy_meas, color="k", label="Measured")
        axC.plot(tms, residual, color="0.5", ls=":",
                 label="Residual")
        axC.axhline(0, color="k", lw=0.6, alpha=0.4)
        axC.set_xlabel("time (ms)")
        axC.set_ylabel("Cytosol cumulative Ca (ions)")
        legend_below(axC, ncol=2)
        axC.margins(x=0)

        print(f"ER budget closure:  max|Σ-measured| = {np.abs(er_sum-er_meas).max():.1f} ions")
        print(f"Cyto residual (end): {residual[-1]:+.0f} ions "
              f"(measured ΔN end = {cy_meas[-1]:+.0f})")

        out = os.path.join(OUTDIR, "budget.png")
        save(figE, panel_path(out, "er"))
        save(figC, panel_path(out, "cyto"))


    plot_fluxes()
    plot_ip3r_comparison()
    plot_budget()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    for name, function in (('traces', run_traces), ('calcium', run_calcium), ('flux', run_flux)):
        p = sub.add_parser(name)
        p.add_argument('model_dir', type=Path)
        p.add_argument('--output-dir', type=Path)
        p.set_defaults(function=function)
        if name == 'traces':
            p.add_argument('--seed', type=int, default=1)
        if name == 'calcium':
            p.add_argument('--ode-csv', default=None, metavar='PATH', help='Optional ODE traces.csv for an ER overlay.')
    args = parser.parse_args()
    apply_style()
    args.function(args)


if __name__ == '__main__':
    main()
