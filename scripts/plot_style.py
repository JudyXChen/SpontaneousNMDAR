"""Project palette and small, single-panel Matplotlib figures."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FIGSIZE = (3, 3)
DPI = 300
LINEWIDTH = 1.6
MARKER_SIZE = 22

# Main and dark tones are for lines; light tones are for fills and bands.
# Each line tone has at least 3:1 contrast against white.
FAMILIES = {
    "blue": {"light": "#B3D2E6", "main": "#4F93BF", "dark": "#2E6F9D",
             "deep": "#174F78"},
    "purple": {"light": "#D7B9E2", "main": "#AF72C3", "dark": "#844A9D"},
    "orange": {"light": "#F2CDA9", "main": "#CC762B", "dark": "#9A4E16"},
}

# Paired categorical colors for NMDAR condition comparisons. Bright colors
# identify the original model (and are reserved for future 3-receptor/cluster
# distance cases); dark colors identify the calibrated model (and the current
# 10-receptor/cluster distance cases).
CONDITION_FAMILIES = {
    "blue": {"bright": "#48CAE4", "dark": "#023047"},
    "purple": {"bright": "#9D4EDD", "dark": "#5A189A"},
    "orange": {"bright": "#FFB703", "dark": "#FB8500"},
}


def shade(family, tone="main"):
    """Return a named project color, e.g. shade('blue', 'dark')."""
    return FAMILIES[family][tone]


def condition_color(family, tone):
    """Return a bright or dark categorical NMDAR condition color."""
    return CONDITION_FAMILIES[family][tone]


def colors(n):
    """Return distinct line colors, using dark tones after the main hues."""
    if not 0 <= n <= 6:
        raise ValueError("n must be between 0 and 6; distinguish more series by line style")
    sequence = [shade(family, tone) for tone in ("main", "dark")
                for family in FAMILIES]
    return sequence[:n]


def blue_colormap():
    """Sequential blue map for ordered heatmap values."""
    return LinearSegmentedColormap.from_list(
        "spont_nmdar_blue", ("#FFFFFF", shade("blue", "light"),
                             shade("blue"), shade("blue", "dark")))


def apply_style():
    plt.switch_backend("Agg")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 1.0,
        "axes.prop_cycle": plt.cycler(color=colors(6)),
        "lines.linewidth": LINEWIDTH,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": DPI,
    })


def figure():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    return fig, ax


def legend_below(ax, ncol=1, **kwargs):
    """Keep legend labels outside the data rectangle on a small figure."""
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                     ncol=ncol, frameon=False, **kwargs)


def condition_output_dir(model_dir):
    """Group individual conditions under their kinetics/protocol branch."""
    from scripts.paths import results_for
    return results_for(model_dir)


def save(fig, path, *, legend_name="auto"):
    """Save a clean panel in fig/ and a viewing copy with a useful legend.

    Label the plotted series to opt into a legend. A lone labeled series needs
    no legend. Pass legend_name=None when a shared legend is exported separately.
    Viewing copies may be taller to keep legends outside a full-size data panel.
    """
    path = Path(path).with_suffix(".png")
    path.parent.mkdir(parents=True, exist_ok=True)
    manuscript = path.parent / "fig" / path.name
    manuscript.parent.mkdir(parents=True, exist_ok=True)
    existing = list(fig.legends) + [ax.get_legend() for ax in fig.axes
                                   if ax.get_legend() is not None]
    handles, labels = [], []
    for ax in fig.axes:
        hs, ls = ax.get_legend_handles_labels()
        for handle, label in zip(hs, ls):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    useful = len(labels) > 1 or bool(fig.legends)
    for legend in existing:
        legend.set_visible(False)
    fig.tight_layout()
    fig.savefig(manuscript, dpi=DPI)
    print(f"Saved -> {manuscript}")
    if useful:
        if existing:
            for legend in existing:
                legend.set_visible(True)
        else:
            legend_below(fig.axes[0], handles=handles, labels=labels)
        fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
        if legend_name is not None:
            if not handles:  # Explicit figure legends, e.g. conductance.
                for legend in fig.legends:
                    handles.extend(legend.legend_handles)
            name = f"{path.stem}_legend.png" if legend_name == "auto" else legend_name
            save_legend(handles, path.parent / name)
    else:
        fig.savefig(path, dpi=DPI)
    print(f"Saved -> {path}")
    plt.close(fig)


def save_legend(handles, path, figsize=(2.8, 0.6), ncol=1, **kwargs):
    """Export a reusable transparent PNG legend in the manuscript fig/ folder."""
    path = Path(path).with_suffix(".png")
    if path.parent.name != "fig":
        path = path.parent / "fig" / path.name
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=figsize)
    fig.legend(handles=handles, loc="center", ncol=ncol,
               frameon=False, **kwargs)
    fig.savefig(path, dpi=DPI, transparent=True,
                bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {path}")
    plt.close(fig)


def panel_path(path, suffix):
    path = Path(path)
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")
