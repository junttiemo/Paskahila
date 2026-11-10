#!/usr/bin/env python3
"""Compare the production seed 1 triple with the Grid B seed 1 triple."""

import argparse
import math
import os
from pathlib import Path
import runpy
import tempfile

import numpy as np


WORKDIR = Path(__file__).resolve().parent
QUEUE = runpy.run_path(str(WORKDIR / "run-production-seed1.py"))
GRIDS = ("production", "Grid B")
LATTICES = (QUEUE["PRODUCTION"], QUEUE["COARSE"])
MODES = (*QUEUE["BASELINE"], "induced")
DT = QUEUE["DT"]
REPORTED = np.arange(111) * 8 / 10


def native_profiles(paths, read_profile):
    fields = np.array([[read_profile(path / f"EMTDynT{t:g}ID1.txt") for t in REPORTED] for path in paths])
    induced = fields[0] - fields[1] - fields[2]
    return np.concatenate((fields, induced[None]), axis=0)


def longitudinal_mapping():
    spacings = [lattice.spacings[2] for lattice in LATTICES]
    z = [np.arange(lattice.dimensions[2]) * dz for lattice, dz in zip(LATTICES, spacings)]
    fine, coarse = zeta = [sites - lattice.crossing[1] for sites, lattice in zip(z, LATTICES)]
    fine_half, coarse_half = (spacing / 2 for spacing in spacings)
    lower = np.maximum(coarse[:, None] - coarse_half, fine - fine_half)
    upper = np.minimum(coarse[:, None] + coarse_half, fine + fine_half)
    overlap = np.maximum(upper - lower, 0)
    covered = np.isclose(overlap.sum(axis=1), spacings[1], rtol=0, atol=1e-12)
    return z, zeta, covered, overlap[covered] / spacings[1]


def time_mapping():
    fine, coarse = [REPORTED + DT / 2 - lattice.crossing[0] for lattice in LATTICES]
    # Evolved times start at the second saved frame. Keep the mixed initial frame separate.
    targets = np.flatnonzero((coarse >= fine[1]) & (coarse <= fine[-1]))
    upper = np.searchsorted(fine, coarse[targets])
    lower = upper - 1
    weight = (coarse[targets] - fine[lower]) / (fine[upper] - fine[lower])
    if len(targets) != 109 or lower[0] < 1 or np.any((weight < 0) | (weight > 1)):
        raise ValueError("Unexpected common evolved time range")
    return coarse, targets, lower, upper, weight


def write_audit(write_csv, output, runs):
    index = []
    for grid, paths in zip(GRIDS, runs):
        for mode, path in zip(MODES, paths):
            gauss = np.loadtxt(path / "GausCheck1.txt")
            peak = gauss[gauss[:, 1].argmax()]
            rms = gauss[gauss[:, 2].argmax()]
            if grid == "production":
                manifest = path / "manifest.json"
                record = QUEUE["read_json"](manifest)
                keys = ("binary_sha256", "wall_seconds", "peak_rss_kib")
            else:
                manifest = path / "MANIFEST.txt"
                record = QUEUE["GRID_B"]["read_fields"](manifest)
                record["executable_hash"] = record["executable_hash"].split()[0]
                keys = ("executable_hash", "elapsed_seconds", "peak_rss_kb")
            common = ("source_revision", "started", "finished", "exit_status", "actual_steps", "actual_endpoint")
            index.append((grid, mode, manifest.relative_to(WORKDIR).as_posix(),
                          *(record[key] for key in (*common, *keys)), peak[1], peak[0], rms[2], rms[0]))
    write_csv(output / "manifest-index.csv", (
        "grid", "mode", "manifest", "source_revision", "started", "finished", "exit_status", "actual_steps",
        "actual_endpoint", "binary_sha256", "wall_seconds", "peak_rss_kib",
        "peak_gauss_max", "peak_gauss_max_time", "peak_gauss_rms", "peak_gauss_rms_time",
    ), index)


def write_metrics(write_csv, output, components, time_indices, s, ref, candidate):
    delta = candidate - ref
    ref_rms = np.sqrt(np.mean(ref**2, axis=2))
    abs_rms = np.sqrt(np.mean(delta**2, axis=2))
    max_abs = np.max(np.abs(delta), axis=2)
    threshold = np.maximum(1e-30, 1e-12 * ref_rms.max(axis=1))
    rows = []
    for m, mode in enumerate(MODES):
        for i, frame in enumerate(time_indices):
            for c, component in enumerate(components):
                denom = ref_rms[m, i, c]
                rms = abs_rms[m, i, c]
                rows.append((
                    "evolved" if frame else "initial_mixed_time", REPORTED[frame], f"{s[frame]:.3f}" if frame else "",
                    mode, component, ref.shape[2], denom, rms, max_abs[m, i, c], rms / denom if denom else "",
                    threshold[m, c], int(denom <= threshold[m, c]),
                ))
    write_csv(output / "profile-differences.csv", (
        "frame_type", "grid_b_t_reported", "s", "mode", "component", "fully_covered_bins", "grid_b_rms",
        "absolute_rms_difference", "maximum_absolute_difference", "relative_rms_difference",
        "near_zero_threshold", "near_zero_grid_b_reference",
    ), rows)


def plot_profiles(plt, output, name, description, labels, x, fields, xlabel):
    fig, axes = plt.subplots(4, 4, figsize=(13, 8.5), sharex=True)
    for m, mode in enumerate(MODES):
        for c, component in enumerate(labels):
            ax = axes[m, c]
            ax.plot(x[0], fields[0][m, :, c], color="#a54b15", linewidth=1.2, label=GRIDS[0])
            ax.plot(x[1], fields[1][m, :, c], color="#26568a", linewidth=1.1, linestyle="--", label=GRIDS[1])
            ax.grid(alpha=0.24)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
            if m == 0:
                ax.set_title(component)
            if c == 0:
                ax.set_ylabel(mode.replace("-", " "))
    fig.suptitle(description, y=0.98)
    fig.supxlabel(xlabel, y=0.065)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=2, frameon=False)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.89, bottom=0.13, wspace=0.37, hspace=0.38)
    fig.savefig(output / f"{name}.png", dpi=145)
    plt.close(fig)


def compare_bulk(plt, write_csv, output, runs):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    rows = []
    read_rows = QUEUE["GRID_B"]["read_numeric_rows"]
    for grid, lattice, paths, color, style in zip(GRIDS, LATTICES, runs, ("#a54b15", "#26568a"), ("-", "--")):
        volume = math.prod(lattice.spacings)
        energies = []
        times = None
        for path in paths:
            initial, *evolved = read_rows(path / "EnergyID1.txt", 4, initial_width=6)
            saved_times = [initial[0], *(row[1] for row in evolved)]
            if times is not None and times != saved_times:
                raise ValueError(f"Bulk time axes differ for {grid} in {path}")
            times = saved_times
            energies.append(volume * np.array([initial[1], *(row[2] for row in evolved)]))
        energies.append(energies[0] - energies[1] - energies[2])
        for ax, mode, values in zip(axes.flat, MODES, energies):
            rows.extend((grid, mode, t, volume, energy) for t, energy in zip(times, values))
            ax.plot(times, values, color=color, linestyle=style, label=grid)
            ax.set_title(mode.replace("-", " "))
            ax.grid(alpha=0.24)
    write_csv(output / "integrated-energy.csv", ("grid", "mode", "t_reported", "cell_volume", "integrated_T00"), rows)
    fig.suptitle("Seed 1 integrated bulk energy")
    fig.supxlabel("Saved bulk time (simulation units)", y=0.07)
    fig.supylabel("Cell volume × bulk site sum")
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=2, frameon=False)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.16, hspace=0.3)
    fig.savefig(output / "integrated-energy.png", dpi=145)
    plt.close(fig)


def compare(args, ensemble):
    write_csv, plt = ensemble["write_csv"], ensemble["plt"]
    runs = [[path.resolve() for path in paths] for paths in (args.production, args.grid_b)]
    for lattice, paths in zip(LATTICES, runs):
        for path in paths:
            QUEUE["validate_info"](path, lattice)
    native = [native_profiles(paths, ensemble["read_emt_profile"]) for paths in runs]
    z, zeta, covered, weights = longitudinal_mapping()
    s, targets, low, high, alpha = time_mapping()
    coarse_zeta = zeta[1][covered]
    b_index = np.flatnonzero(covered)
    b_bins, p_cells = np.nonzero(weights)
    write_csv(args.output / "z-mapping.csv", (
        "grid_b_z_index", "grid_b_z", "grid_b_zeta", "production_z_index", "production_z", "production_zeta",
        "production_weight_in_bin",
    ), ((b_index[b], z[1][b_index[b]], coarse_zeta[b], p, z[0][p], zeta[0][p], weights[b, p])
        for b, p in zip(b_bins, p_cells)))
    write_csv(args.output / "time-mapping.csv", (
        "frame_type", "grid_b_t_reported", "s", "production_t_reported_lower", "production_t_reported_upper",
        "upper_time_weight",
    ), [("initial_mixed_time", 0, "", 0, 0, ""),
        *(("evolved", REPORTED[b], f"{s[b]:.3f}", REPORTED[l], REPORTED[h], f"{a:.8g}")
          for b, l, h, a in zip(targets, low, high, alpha))])

    rebinned = np.einsum("bp,mtpq->mtbq", weights, native[0], optimize=True)
    interpolated = rebinned[:, low] + (rebinned[:, high] - rebinned[:, low]) * alpha[None, :, None, None]
    production = np.concatenate((rebinned[:, :1], interpolated), axis=1)
    time_indices = np.r_[0, targets]
    baseline = native[1][:, time_indices][:, :, covered]
    write_metrics(write_csv, args.output, ensemble["COMPONENTS"], time_indices, s, baseline, production)
    write_audit(write_csv, args.output, runs)
    compare_bulk(plt, write_csv, args.output, runs)

    after = [int(np.searchsorted(REPORTED + DT / 2, lattice.crossing[0])) for lattice in LATTICES]
    snapshots = (("initial", 0, 0), ("before-crossing", after[0] - 1, after[1] - 1),
                 ("after-crossing", *after), ("t24", 30, 30), ("t40", 50, 50), ("t88", 110, 110))
    for name, p, b in snapshots:
        index = int(np.searchsorted(time_indices, b))
        plot_profiles(plt, args.output, f"native-{name}",
                      f"Native seed 1 profiles | Production t={REPORTED[p]:g}, Grid B t={REPORTED[b]:g}",
                      ensemble["LABELS"], z, [native[0][:, p], native[1][:, b]], "Saved z (simulation units)")
        phase = f"s={s[b]:.3f} at Grid B t={REPORTED[b]:g}" if b else "Initial mixed-time frame"
        plot_profiles(plt, args.output, f"aligned-{name}", f"Aligned seed 1 profiles | {phase}", ensemble["LABELS"],
                      (coarse_zeta, coarse_zeta), (production[:, index], baseline[:, index]),
                      "Position relative to nominal crossing ζ (simulation units)")

    QUEUE["write_json"](args.output / "comparison.json", {
        "seed": 1, "dT": DT, "modes": MODES, "components": ensemble["COMPONENTS"],
        "grid": {grid: {"crossing_t_z": lattice.crossing, "dimensions": lattice.dimensions,
                        "spacings": lattice.spacings} for grid, lattice in zip(GRIDS, LATTICES)},
        "common_z_bins": int(covered.sum()), "matched_evolved_times": len(targets),
        "evolved_time_coordinate": "s = t_reported + dT/2 - t_crossing",
        "space_coordinate": "zeta = z_index * az - z_crossing; overlap weights use cell edges",
        "near_zero_rule": "RMS(Grid B) <= max(1e-30, 1e-12 * peak RMS for this mode and component)",
        "interpretation": "The comparison combines mesh and transverse charge realization effects for these events.",
    })
    print(f"Compared {len(targets)} evolved times and the initial frame in {covered.sum()} Grid B bins", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-imports", action="store_true")
    parser.add_argument("--production", type=Path, nargs=3, metavar="RUN")
    parser.add_argument("--grid-b", type=Path, nargs=3, metavar="RUN")
    parser.add_argument("--output", type=Path, help="existing comparison output directory")
    args = parser.parse_args()
    if not args.check_imports and not (args.production and args.grid_b and args.output and args.output.is_dir()):
        parser.error("Provide both triples in collision, left-only, right-only order and an existing --output directory")
    with tempfile.TemporaryDirectory(prefix="production-mpl-") as cache:
        os.environ["MPLCONFIGDIR"] = cache
        ensemble = runpy.run_path(str(WORKDIR / "analyze-grid-b-ensemble.py"))
        if args.check_imports:
            print(f"Analysis imports OK | NumPy {np.__version__} | Matplotlib {ensemble['matplotlib'].__version__}",
                  flush=True)
        else:
            compare(args, ensemble)


if __name__ == "__main__":
    main()
