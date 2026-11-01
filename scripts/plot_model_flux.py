"""Ca2+ fluxes and compartment budgets for the IP3R_homog run (5 seeds).

All fluxes are reconstructed from the recorded molecule-count observables
(react_data/seed_*/*.dat) and the model's rate laws, in a common unit of
ions per second (ions s^-1), using the **cytosol-positive** sign convention:

    flux INTO the cytosol  > 0        flux OUT of the cytosol  < 0
      ER  -> cytosol  (+)               cytosol -> ER   (-)
      ECS -> cytosol  (+)               cytosol -> ECS  (-)

Separate 3 x 3 inch figures are written to the matching model branch under
results/:
  flux_<pathway>.png -- one signed flux per file (mean +/-1 SD)
  ip3r_flux.png               -- direct and reconstructed IP3R flux
  budget_{er,cyto}.png -- one compartment per file.

Flux definitions (see model.bngl / subsystem.py):
    J_VSCC    = vscc_rate(t) * N_VSCC                     ECS -> cyto  (+)
    J_NMDAR   = nmdar_rate(t) * N_NMDAR                   ECS -> cyto  (+)
    J_PMCA,L  = k_PMCA_leak * (N_PMCA_tot - N_PMCA1)      ECS -> cyto  (+)
    J_NCX,L   = k_NCX_leak  * (N_NCX_tot  - N_NCX1)       ECS -> cyto  (+)
    J_PMCA    = k_PMCA_3 * N_PMCA1                        cyto -> ECS  (-)
    J_NCX     = k_NCX_3  * N_NCX1                         cyto -> ECS  (-)
    J_SERCA,L = k_SERCA_leak * N_CaER                     ER  -> cyto  (+)
    J_RyR     = k_FluxRyR      * N_RyR_O   *([CaER]-[Cacyto])  ER -> cyto (+)
    J_IP3R    = k_IP3R_Ca_flux * N_IP3R_110 *[CaER]       ER  -> cyto  (+)
    J_SERCA   = -(d N_CaER/dt) + J_SERCA,L + J_RyR + J_IP3R    cyto -> ER (-)
                (uptake, by ER mass balance -- SERCA cycle states were not recorded)

Concentrations use the volumes the model implies: V = N(t=0)/(N_A * IC_conc).
Pump totals use placement density x membrane area (see constants).
"""
import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for
from scripts.plot_style import condition_output_dir, apply_style, figure, legend_below, panel_path, save, shade
apply_style()

# --------------------------------------------------------------------------
# constants (model.bngl parameters + instantiation.py placement)
# --------------------------------------------------------------------------
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

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument('model_dir', type=Path, help='Model leaf; recorded data is resolved under results/')
ap.add_argument('--output-dir', type=Path, default=None)
ARGS = ap.parse_args()
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


if __name__ == "__main__":
    plot_fluxes()
    plot_ip3r_comparison()
    plot_budget()
