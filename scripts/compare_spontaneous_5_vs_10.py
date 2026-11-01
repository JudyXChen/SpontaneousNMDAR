"""Plot and validate the matched spontaneous 5- versus 10-NMDAR seed-1 run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import sys

import numpy as np
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import results_for

from scripts.plot_style import (  # noqa: E402
    apply_style,
    condition_color,
    figure,
    save,
    save_legend,
)


RELEASE_TIME_S = 0.010
EXPECTED_END_TIME_S = 0.045
EXPECTED_RELEASE_ITERATION = 100_000
TIME_TOLERANCE_S = 1e-10

CONDITIONS = {
    "10_NMDAR": {
        "label": "10 NMDAR",
        "color": condition_color("blue", "dark"),
        "leaf": ROOT / "model/calcium_total/spont/1_cluster_10_NMDAR",
    },
    "5_NMDAR": {
        "label": "5 NMDAR",
        "color": condition_color("orange", "dark"),
        "leaf": ROOT / "model/calcium_total/spont/1_cluster_5_NMDAR",
    },
}

TRACE_FILES = {
    "nmdar_ca_entry_firings": "NMDAR_Ca_flux_World.dat",
    "open_nmdar_count": "NMDAR_World.dat",
    "spine_cytosolic_free_ca_ions": "Ca2_Mushroom_transparent.dat",
}

PANELS = {
    "NMDAR_Ca_entry": (
        "nmdar_ca_entry_firings",
        "Cumulative modeled NMDAR\n" + r"Ca$^{2+}$ entry (firings)",
    ),
    "open_NMDAR": ("open_nmdar_count", "Open NMDAR count"),
    "spine_cytosolic_Ca_raw": (
        "spine_cytosolic_free_ca_ions",
        r"Spine cytosolic free Ca$^{2+}$ (ions)",
    ),
    "spine_cytosolic_Ca_delta": (
        "spine_cytosolic_free_ca_delta_ions",
        r"Baseline-subtracted spine Ca$^{2+}$ (ions)",
    ),
}


def load_trace(path: Path) -> np.ndarray:
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(f"Expected a two-column trace: {path}")
    if not np.isclose(values[0, 0], 0, atol=TIME_TOLERANCE_S):
        raise ValueError(f"Trace does not start at 0 s: {path}")
    if not np.isclose(
        values[-1, 0], EXPECTED_END_TIME_S, atol=TIME_TOLERANCE_S
    ):
        raise ValueError(
            f"Incomplete trace ends at {values[-1, 0]:g} s, "
            f"expected {EXPECTED_END_TIME_S:g} s: {path}"
        )
    if np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError(f"Trace time must be strictly increasing: {path}")
    return values


def reaction_probability_warnings(leaf: Path) -> list[str]:
    report_path = leaf / "reports/warnings_report_00001.txt"
    if report_path.exists():
        report_text = report_path.read_text(errors="replace").strip()
        return [
            block.strip()
            for block in re.split(r"(?=^Warning:)", report_text, flags=re.MULTILINE)
            if block.strip()
        ]

    warnings: list[str] = []
    for path in (leaf / "logs/seed_1.log",):
        if not path.exists():
            continue
        for line in path.read_text(errors="replace").splitlines():
            normalized = line.strip()
            lowered = normalized.lower()
            if lowered.startswith("warning") and "probability" in lowered:
                if normalized not in warnings:
                    warnings.append(normalized)
    return warnings


def validate_run(leaf: Path) -> dict[str, object]:
    report_path = leaf / "reports/run_report_00001.txt"
    log_path = leaf / "logs/seed_1.log"
    report_text = report_path.read_text(errors="replace")
    log_text = log_path.read_text(errors="replace")

    releases = [
        (match.group("description"), int(match.group("iteration")))
        for match in re.finditer(
            r'^Released (?P<description>.+?) at iteration (?P<iteration>\d+)\.$',
            log_text,
            flags=re.MULTILINE,
        )
    ]
    glutamate_iterations = [
        iteration
        for description, iteration in releases
        if '500 glu from "Release_Site_11"' in description
    ]
    other_release_iterations = [
        iteration
        for description, iteration in releases
        if 'glu from "Release_Site_11"' not in description
    ]

    finished = report_text.rstrip().endswith("FINISHED")
    glutamate_at_expected_iteration = (
        glutamate_iterations == [EXPECTED_RELEASE_ITERATION]
    )
    other_species_at_zero = bool(other_release_iterations) and all(
        iteration == 0 for iteration in other_release_iterations
    )
    if not finished:
        raise ValueError(f"Run report does not end in FINISHED: {report_path}")
    if not glutamate_at_expected_iteration:
        raise ValueError(
            f"Expected one glutamate release at iteration "
            f"{EXPECTED_RELEASE_ITERATION}, found {glutamate_iterations}: {log_path}"
        )
    if not other_species_at_zero:
        raise ValueError(f"Non-glutamate initialization was not at iteration 0: {log_path}")

    return {
        "run_report_ends_finished": finished,
        "glutamate_release_iteration": glutamate_iterations[0],
        "glutamate_release_at_expected_iteration": glutamate_at_expected_iteration,
        "all_other_logged_releases_at_iteration_zero": other_species_at_zero,
        "logged_release_count": len(releases),
    }


def save_panel(fig, output_stem: Path) -> None:
    save(fig, output_stem.with_suffix(".png"), legend_name=None)


def scalar(value: float) -> int | float:
    value = float(value)
    if np.isclose(value, round(value), atol=1e-12):
        return int(round(value))
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            ROOT
            / "results/calcium_total/spont/"
            "compare_1cluster"
        ),
    )
    parser.add_argument("--raw-only", action="store_true",
                        help="Omit the optional baseline-subtracted calcium panel.")
    args = parser.parse_args()

    apply_style()
    traces: dict[tuple[str, str], np.ndarray] = {}
    baseline_by_condition: dict[str, float] = {}
    reference_time: np.ndarray | None = None
    validation: dict[str, dict[str, object]] = {}

    for condition, config in CONDITIONS.items():
        leaf = results_for(config["leaf"])
        validation[condition] = validate_run(leaf)
        for quantity, file_name in TRACE_FILES.items():
            values = load_trace(leaf / "react_data/seed_00001" / file_name)
            if reference_time is None:
                reference_time = values[:, 0]
            elif not np.array_equal(values[:, 0], reference_time):
                raise ValueError(
                    f"Time grid differs from the other traces: {leaf / file_name}"
                )
            traces[(condition, quantity)] = values

    assert reference_time is not None
    shifted_time_s = reference_time - RELEASE_TIME_S
    baseline_mask = (
        (shifted_time_s >= -0.005 - TIME_TOLERANCE_S)
        & (shifted_time_s < -TIME_TOLERANCE_S)
    )
    post_release_mask = shifted_time_s >= -TIME_TOLERANCE_S
    if not np.any(baseline_mask):
        raise ValueError("The -5 to 0 ms baseline window contains no samples")

    for condition in CONDITIONS:
        calcium = traces[(condition, "spine_cytosolic_free_ca_ions")][:, 1]
        baseline = float(np.mean(calcium[baseline_mask]))
        baseline_by_condition[condition] = baseline
        traces[(condition, "spine_cytosolic_free_ca_delta_ions")] = np.column_stack(
            (reference_time, calcium - baseline)
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "traces.csv"
    csv_columns = [
        "nmdar_ca_entry_firings",
        "open_nmdar_count",
        "spine_cytosolic_free_ca_ions",
        "spine_cytosolic_free_ca_delta_ions",
    ]
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["simulation_time_s", "time_from_glutamate_release_ms"]
            + [
                f"{quantity}_{condition}"
                for quantity in csv_columns
                for condition in CONDITIONS
            ]
        )
        for index, simulation_time_s in enumerate(reference_time):
            writer.writerow(
                [
                    f"{simulation_time_s:.8g}",
                    f"{shifted_time_s[index] * 1e3:.8g}",
                ]
                + [
                    f"{traces[(condition, quantity)][index, 1]:.8g}"
                    for quantity in csv_columns
                    for condition in CONDITIONS
                ]
            )
    print(f"Saved -> {csv_path}")

    for panel_name, (quantity, ylabel) in PANELS.items():
        if args.raw_only and panel_name == "spine_cytosolic_Ca_delta":
            continue
        fig, ax = figure()
        for condition, config in CONDITIONS.items():
            values = traces[(condition, quantity)]
            ax.plot(
                shifted_time_s * 1e3,
                values[:, 1],
                color=config["color"],
                drawstyle="steps-post",
                label=config["label"],
            )
        ax.set(
            xlabel="Time from glutamate release (ms)",
            ylabel=ylabel,
            xlim=(-10, 35),
        )
        ax.set_xticks((-10, 0, 10, 20, 30))
        if panel_name in {"NMDAR_Ca_entry", "open_NMDAR"}:
            ax.set_ylim(bottom=-0.05)
        save_panel(fig, args.output_dir / panel_name)

    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [],
            [],
            color=config["color"],
            linewidth=1.6,
            label=config["label"],
        )
        for config in CONDITIONS.values()
    ]
    save_legend(
        handles,
        args.output_dir / "legend.png",
        figsize=(2.8, 0.55),
        ncol=2,
    )

    metrics: dict[str, object] = {
        "protocol": {
            "seed": 1,
            "time_step_s": 1e-7,
            "iterations": 450_000,
            "glutamate_release_time_s": RELEASE_TIME_S,
            "post_release_duration_s": 0.035,
            "baseline_window_ms_from_release": [-5, 0],
            "baseline_window_excludes_release_sample": True,
            "integrated_open_occupancy_method": (
                "trapezoidal integral of sampled open count from 0 to 35 ms"
            ),
        },
        "validation": validation,
        "conditions": {},
    }
    sample_interval_s = round(float(np.median(np.diff(reference_time))), 12)
    for condition, config in CONDITIONS.items():
        open_count = traces[(condition, "open_nmdar_count")][:, 1]
        ca_entry = traces[(condition, "nmdar_ca_entry_firings")][:, 1]
        calcium = traces[(condition, "spine_cytosolic_free_ca_ions")][:, 1]
        calcium_delta = traces[
            (condition, "spine_cytosolic_free_ca_delta_ions")
        ][:, 1]
        warnings = reaction_probability_warnings(config["leaf"])
        integrated_open_occupancy_s = float(
            (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(
                open_count[post_release_mask],
                reference_time[post_release_mask],
            )
        )
        metrics["conditions"][condition] = {
            "label": config["label"],
            "trace_sample_count": int(reference_time.size),
            "trace_sample_interval_s": sample_interval_s,
            "trace_end_time_s": float(reference_time[-1]),
            "peak_open_nmdar_count": scalar(np.max(open_count[post_release_mask])),
            "integrated_sampled_open_occupancy_receptor_s": (
                integrated_open_occupancy_s
            ),
            "integrated_sampled_open_occupancy_receptor_ms": (
                integrated_open_occupancy_s * 1e3
            ),
            "final_cumulative_nmdar_ca_entry_firings": scalar(ca_entry[-1]),
            "baseline_spine_cytosolic_free_ca_ions": float(
                baseline_by_condition[condition]
            ),
            "peak_raw_spine_cytosolic_free_ca_ions": scalar(
                np.max(calcium[post_release_mask])
            ),
            "peak_baseline_subtracted_spine_ca_ions": float(
                np.max(calcium_delta[post_release_mask])
            ),
            "mcell_reaction_probability_warning_count": len(warnings),
            "mcell_reaction_probability_warnings": warnings,
        }

    metrics_path = args.output_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"Saved -> {metrics_path}")


if __name__ == "__main__":
    main()
