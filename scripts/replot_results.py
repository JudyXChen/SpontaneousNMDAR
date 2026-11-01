"""Regenerate the existing result figures using the shared PNG export layout.

Run with the Blender Python or another Python with NumPy, Matplotlib and SciPy.
No simulations are run. Original and calibrated conductance curves are read
from their respective saved CSVs. Other panels use recorded model outputs.
"""

from pathlib import Path
import os
import argparse
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.paths import ARCHIVE, model_family, result_family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-archive', action='store_true',
                        help='also regenerate old-kinetics figures and comparisons')
    args = parser.parse_args()
    families = ['calcium_total'] + (['NMDAR_old'] if args.include_archive else [])
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/spont-nmdar-mpl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/spont-nmdar-cache")
    os.environ.setdefault("PYTHONPYCACHEPREFIX", "/tmp/spont-nmdar-pycache")
    results = ROOT / "results"
    # Only move figure-condition folders; isolated simulation runs stay in place.
    for kinetics in families:
        branch = result_family(kinetics)
        for protocol in ("evoked", "spont", "mixed"):
            parent = branch / protocol
            if not parent.exists():
                continue
            for leaf in sorted(parent.iterdir()):
                if (leaf.is_dir() and not leaf.name.startswith("compare_")
                        and leaf.name != "conditions"
                        and (model_family(kinetics) / protocol / leaf.name).is_dir()
                        and any(leaf.rglob("*.png"))):
                    target = parent / "conditions" / leaf.name
                    if target.exists():
                        raise FileExistsError(f"Both old and grouped condition folders exist: {leaf}")
                    target.parent.mkdir(exist_ok=True)
                    shutil.move(str(leaf), target)
                    print(f"Moved -> {target}", flush=True)

    # Track the exact existing panels and require freshly generated replacements
    # before deleting any superseded image files.
    figure_roots = [result_family(family) for family in families]
    if args.include_archive:
        figure_roots.append(ARCHIVE / 'comparison')
    panels = [p for folder in figure_roots for p in folder.rglob('*.png')
              if p.parent.name != 'fig' and 'legend' not in p.stem]
    started = time.time()
    tasks = []
    for kinetics in families:
        output = result_family(kinetics) / "conductance"
        if (output / "comparison.csv").exists():
            tasks.append(("nmdar/compare_conductance.py", "--from-csv",
                          "--output-dir", output))
    old = ARCHIVE / "results/evoked"
    if args.include_archive and (old / "conditions/1_cluster_10_NMDAR").exists():
        tasks.extend((name,) for name in (
            "plot_cluster_validation.py", "plot_cluster_layout.py",
            "plot_calcium_positions.py", "plot_calcium_depth.py"))
    for kinetics in families:
        for protocol in ("evoked", "spont", "mixed"):
            grouped = result_family(kinetics) / protocol / "conditions"
            if not grouped.exists():
                continue
            for output in sorted(grouped.iterdir()):
                if not output.is_dir():
                    continue
                model = model_family(kinetics) / protocol / output.name
                if (output / "ca_er.png").exists():
                    tasks.append(("plot_model_calcium.py", model))
                if (output / "flux_nmdar.png").exists():
                    tasks.append(("plot_model_flux.py", model))
                if (kinetics == "calcium_total" and protocol == "evoked"
                        and (output / "psd_layout.png").exists()):
                    tasks.append(("plot_evoked_cluster_case.py", model))
    calibrated = results / "calcium_total/evoked"
    comparisons = (
        (calibrated / "compare_1or2cluster", "compare_cluster_calibration.py",
         ("--calibrated-only",)),
        (calibrated / "compare_2cluster_dist", "plot_distance.py", ()),
        (calibrated / "compare_2cluster3NMDAR_dist",
         "plot_distance.py", ("--receptors", "3")),
        (ARCHIVE / "comparison/evoked/"
         "compare_1_cluster_10_NMDAR_vs_2_cluster_100nm_seed1",
         "compare_cluster_calibration.py", ()),
        (results / "calcium_total/spont/compare_1cluster",
         "compare_spontaneous_5_vs_10.py",
         (() if (results / "calcium_total/spont/compare_1cluster/"
                 "spine_cytosolic_Ca_delta.png").exists() else ("--raw-only",))),
    )
    for output, script, arguments in comparisons:
        if output.exists() and (args.include_archive or not output.is_relative_to(ARCHIVE)):
            tasks.append((script, *arguments, "--output-dir", output))
    fixed_comparisons = (
        (calibrated / "compare_smallcluster", "plot_smallcluster_opening_comparison.py"),
        (ARCHIVE / "comparison/evoked/"
         "compare_1_cluster_10_NMDAR_seed1", "compare_one_cluster_opening.py"),
        (ARCHIVE / "comparison/kinetics/100uM_Glu_1ms",
         "nmdar/analyze_nmdar_residual_current.py"),
    )
    for output, script in fixed_comparisons:
        if output.exists() and (args.include_archive or not output.is_relative_to(ARCHIVE)):
            tasks.append((script,))
    for index, (script, *arguments) in enumerate(tasks, 1):
        print(f"[{index}/{len(tasks)}] {script}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / script),
                        *map(str, arguments)], cwd=ROOT, check=True)

    for panel in panels:
        target = panel.parent.parent / panel.name if panel.parent.name in (
            "traces", "positions") else panel
        clean = target.parent / "fig" / target.name
        if not (target.exists() and clean.exists()
                and target.stat().st_mtime >= started
                and clean.stat().st_mtime >= started):
            raise RuntimeError(f"Missing freshly generated PNG pair for {panel}")
    removed = 0
    for path in list(results.rglob("*")):
        if path.suffix in (".pdf", ".svg"):
            path.unlink()
            removed += 1
        elif path.suffix == ".png" and path.parent.name != "fig":
            if "legend" in path.stem and (path.parent / "fig" / path.name).exists():
                path.unlink()
            elif path.parent.name in ("traces", "positions"):
                if (path.parent.parent / "fig" / path.name).exists():
                    path.unlink()
    for folder in sorted(results.rglob("*"), reverse=True):
        if folder.is_dir() and folder.name in ("traces", "positions"):
            # Finder metadata is not a result artifact.
            metadata = folder / ".DS_Store"
            if metadata.exists():
                metadata.unlink()
            if not any(folder.iterdir()):
                folder.rmdir()
    print(f"Regenerated {len(panels)} existing panels; removed {removed} PDF/SVG exports.")


if __name__ == "__main__":
    main()
