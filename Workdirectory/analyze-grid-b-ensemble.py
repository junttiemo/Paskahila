#!/usr/bin/env python3
"""Audit and analyze the completed eight-seed Grid B EMT ensemble."""

import argparse
import csv
import json
from pathlib import Path
import re
import runpy

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parent
SEEDS = tuple(range(1, 9))
STEPS = ("0.08", "0.008")
MODES = ("collision", "left-only", "right-only", "induced")
COMPONENTS = ("T00", "T0Z", "P_T", "P_L")
LABELS = (r"$T_{00}$", r"$T_{0z}$", r"$P_T$", r"$P_L$")
REPORTED = np.arange(111) * 8 / 10
Z_INDEX = np.arange(512)
Z = Z_INDEX * 0.5
ZETA = Z - 127.25
T_COLL = 19.25
SNAPSHOTS = (24, 30, 50, 110)  # Frames nearest impact and at reported t=24, 40, 88.
FIELD_NAMES = tuple(f"{mode.replace('-', '_')}_{field}" for mode in MODES for field in COMPONENTS)
TIME_COLUMNS = (
    "t_fine_reported", "t_coarse_lower_reported", "t_coarse_upper_reported",
    "t_physical", "s", "coarse_upper_weight",
)
SPACE_COLUMNS = ("z_index", "z", "zeta")
SPACE = np.column_stack((Z_INDEX, Z, ZETA))
SPACE_FORMAT = ("%d", "%.1f", "%.2f")
TIME_FORMAT = ("%.1f", "%.1f", "%.1f", "%.3f", "%.3f", "%.6f")
VALUE_FORMAT = "%.16g"


def write_csv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def read_emt_profile(path):
    t00, txx, tyy, tzz, t0z = np.loadtxt(path, usecols=(1, 2, 3, 4, 7), unpack=True)
    return np.column_stack((t00, t0z, (txx + tyy) / 2, tzz))


def write_profile_csv(path, columns, coordinate_format, chunks):
    fmt = coordinate_format + (VALUE_FORMAT,) * (len(columns) - len(coordinate_format))
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write(",".join(columns) + "\n")
        for times, space, fields in chunks:
            coordinates = np.column_stack((np.repeat(times, len(space), axis=0), np.tile(space, (len(times), 1))))
            rows = np.column_stack((coordinates, fields.reshape(len(coordinates), -1)))
            np.savetxt(stream, rows, fmt=fmt, delimiter=",")
    print(f"Wrote {path.name}", flush=True)


def load_campaign(queue, completed, output):
    # Axes are timestep, seed, run type, frame, z, and component.
    raw = np.empty((2, len(SEEDS), len(MODES), len(REPORTED), len(Z), len(COMPONENTS)))
    index = []
    audit = []
    common_hashes = None
    hash_pattern = re.compile(r"^([0-9a-f]{64})  (SourceCode/\S+)$", re.MULTILINE)

    for dt, seed, mode in queue["JOBS"]:
        run = completed[dt, seed, mode]
        path = run / "MANIFEST.txt"
        manifest = queue["read_fields"](path)
        launch, _, rest = path.read_text(encoding="utf-8").partition("\nbuild_record:\n")
        build, _, _ = rest.partition("\nphysics_and_grid:\n")
        launch_revision = re.search(r"^source_revision=([0-9a-f]+)$", launch, re.MULTILINE)
        source_hashes = {name: digest for digest, name in hash_pattern.findall(launch)}
        build_hashes = {name: digest for digest, name in hash_pattern.findall(build)}
        if not launch_revision or set(source_hashes) != set(queue["BUILD_SOURCES"]) or source_hashes != build_hashes:
            raise ValueError(f"Missing or inconsistent source provenance in {path}")
        if common_hashes is None:
            common_hashes = source_hashes
        if source_hashes != common_hashes:
            raise ValueError(f"Source hashes differ from seed 1 in {path}")

        fields = raw[STEPS.index(dt), int(seed) - 1, MODES.index(mode)]
        for frame, time in enumerate(REPORTED):
            profile = run / f"EMTDynT{time:g}ID{seed}.txt"
            fields[frame] = read_emt_profile(profile)
            physical = f"{time + float(dt) / 2:.3f}" if frame else ""
            if frame == 0:
                alignment = "initial_mixed_time"
            elif dt == STEPS[0]:
                alignment = "coarse_bracket"
            elif frame == 1:
                alignment = "before_coarse_bracket"
            else:
                alignment = "fine_target"
            audit.append((seed, dt, mode, frame, frame * int(manifest["emt_profile_cadence_steps"]),
                          f"{time:.1f}", physical, alignment, len(Z), 14, profile.relative_to(WORKDIR).as_posix()))

        gauss = np.loadtxt(run / f"GausCheck{seed}.txt")
        with (run / f"EnergyID{seed}.txt").open(encoding="utf-8") as stream:
            bulk_t00 = float(stream.readline().split()[1])
        index.append((seed, dt, mode, *manifest["charge_density_multipliers_left_right"].split(),
                      path.relative_to(WORKDIR).as_posix(), manifest["started"], manifest["finished"],
                      launch_revision[1], manifest["source_revision"], manifest["executable_hash"].split()[0],
                      json.dumps(source_hashes, sort_keys=True), manifest["threads"], manifest["actual_steps"],
                      manifest["actual_endpoint"], manifest["elapsed_seconds"], manifest["peak_rss_kb"],
                      manifest["exit_status"], bulk_t00, np.max(np.abs(fields[0, :, 0])),
                      gauss[:, 1].max(), gauss[:, 2].max()))

    write_csv(output / "manifest-index.csv", (
        "seed", "dT", "run_type", "nucleus_left", "nucleus_right", "manifest", "started", "finished",
        "launch_revision", "build_revision", "binary_sha256", "source_sha256", "threads", "actual_steps",
        "actual_endpoint", "elapsed_seconds", "peak_rss_kib", "exit_status", "initial_bulk_T00_site_sum",
        "initial_peak_abs_T00", "peak_gauss_max", "peak_gauss_rms",
    ), index)
    write_csv(output / "profile-time-audit.csv", (
        "seed", "dT", "run_type", "frame", "step", "t_reported", "t_physical", "alignment_role",
        "z_rows", "columns", "profile",
    ), audit)
    raw[:, :, 3] = raw[:, :, 0] - raw[:, :, 1] - raw[:, :, 2]
    return raw


def align_times(raw):
    # Fine frame i lies between coarse frames i-1 and i for i >= 2.
    lower, upper, fine_reported = REPORTED[1:-1], REPORTED[2:], REPORTED[2:]
    physical = np.round(fine_reported + float(STEPS[1]) / 2, 3)
    coarse_offset = float(STEPS[0]) / 2
    weight = (physical - (lower + coarse_offset)) / (upper - lower)
    if not np.all((0 < weight) & (weight < 1)):
        raise ValueError("Fine profile time outside its evolved coarse bracket")
    time_map = np.column_stack((fine_reported, lower, upper, physical, np.round(physical - T_COLL, 3), weight))
    coarse_low, coarse_high = raw[0, :, :, 1:-1], raw[0, :, :, 2:]
    coarse = coarse_low + (coarse_high - coarse_low) * weight[None, None, :, None, None]
    fine = raw[1, :, :, 2:]
    return time_map, coarse, fine


def seed_statistics(fields, axis):
    sd = fields.std(axis=axis, ddof=1)
    return np.stack((fields.mean(axis=axis), sd, sd / np.sqrt(len(SEEDS))), axis=-1)


def write_profiles(output, raw, time_map, coarse, fine, delta):
    write_profile_csv(output / "paired-raw.csv", ("seed", "dT", "t_reported", *SPACE_COLUMNS[:2], *FIELD_NAMES),
                      ("%d", "%.3g", "%.1f", *SPACE_FORMAT[:2]), (
        (np.column_stack((np.full(111, seed), np.full(111, float(dt)), REPORTED)), SPACE[:, :2],
         raw[step, i].transpose(1, 2, 0, 3))
        for step, dt in enumerate(STEPS) for i, seed in enumerate(SEEDS)
    ))

    comparison_fields = tuple(f"{name}_{term}" for name in FIELD_NAMES for term in ("coarse", "fine", "fine_minus_coarse"))
    write_profile_csv(output / "paired-aligned.csv", ("seed", *TIME_COLUMNS, *SPACE_COLUMNS, *comparison_fields),
                      ("%d", *TIME_FORMAT, *SPACE_FORMAT), (
        (np.column_stack((np.full(len(time_map), seed), time_map)), SPACE,
         np.stack((coarse[i], fine[i], delta[i]), axis=-1).transpose(1, 2, 0, 3, 4))
        for i, seed in enumerate(SEEDS)
    ))
    write_profile_csv(output / "initial-frame-comparison.csv", ("seed", "t_reported", *SPACE_COLUMNS, *comparison_fields),
                      ("%d", "%.1f", *SPACE_FORMAT), (
        (np.array([[seed, 0]]), SPACE, np.stack((raw[0, i, :, :1], raw[1, i, :, :1],
                                                raw[1, i, :, :1] - raw[0, i, :, :1]), axis=-1).transpose(1, 2, 0, 3, 4))
        for i, seed in enumerate(SEEDS)
    ))

    statistics = ("mean", "sd", "sem")
    induced_stats = seed_statistics(raw[:, :, 3], axis=1)
    write_profile_csv(output / "ensemble-induced.csv", ("dT", "t_reported", "n_seeds", *SPACE_COLUMNS, *(
        f"{field}_{stat}" for field in COMPONENTS for stat in statistics
    )), ("%.3g", "%.1f", "%d", *SPACE_FORMAT), (
        (np.column_stack((np.full(111, float(dt)), REPORTED, np.full(111, len(SEEDS)))), SPACE, induced_stats[step])
        for step, dt in enumerate(STEPS)
    ))
    difference_stats = seed_statistics(delta, axis=0).transpose(1, 2, 0, 3, 4)
    write_profile_csv(output / "ensemble-timestep-difference.csv", (*TIME_COLUMNS, "n_seeds", *SPACE_COLUMNS, *(
        f"{name}_{stat}" for name in FIELD_NAMES for stat in statistics
    )), (*TIME_FORMAT, "%d", *SPACE_FORMAT), (
        (np.column_stack((time_map, np.full(len(time_map), len(SEEDS)))), SPACE, difference_stats),
    ))


def rms(values, axis=None):
    return np.sqrt(np.mean(values ** 2, axis=axis))


def write_timestep_metrics(output, raw, time_map, coarse, fine):
    cases = [("matched_window", "", "", coarse, fine), ("initial", "0", "", raw[0, :, :, :1], raw[1, :, :, :1])]
    for frame in SNAPSHOTS:
        i = frame - 2
        cases.append(("snapshot", f"{REPORTED[frame]:g}", f"{time_map[i, 3]:.3f}", coarse[:, :, i:i + 1], fine[:, :, i:i + 1]))

    seeds = []
    summary = []
    for window, reported, physical, low, high in cases:
        difference = high - low
        low_rms, high_rms, diff_rms = (rms(values, axis=(2, 3)) for values in (low, high, difference))
        abs_max = np.max(np.abs(difference), axis=(2, 3))
        spread = rms(high.std(axis=0, ddof=1), axis=(1, 2))
        for m, mode in enumerate(MODES):
            for c, component in enumerate(COMPONENTS):
                base, ref, change = low_rms[:, m, c], high_rms[:, m, c], diff_rms[:, m, c]
                prefix = (window, reported, physical, mode, component)
                for i, seed in enumerate(SEEDS):
                    seeds.append((*prefix, seed, base[i], ref[i], change[i], change[i] / base[i] if base[i] else "",
                                  abs_max[i, m, c]))
                pooled_base, pooled_ref, pooled_change = rms(base), rms(ref), rms(change)
                relative = change[base > 0] / base[base > 0]
                summary.append((*prefix, len(SEEDS), pooled_base, pooled_ref, pooled_change,
                                pooled_change / pooled_base if pooled_base else "", spread[m, c],
                                relative.min() if relative.size else "", relative.max() if relative.size else ""))

    prefixes = ("window", "t_fine_reported", "t_physical", "run_type", "component")
    write_csv(output / "seed-timestep-metrics.csv", (*prefixes, "seed", "rms_coarse", "rms_fine", "rms_difference",
              "relative_rms_difference", "max_abs_difference"), seeds)
    write_csv(output / "timestep-summary.csv", (*prefixes, "n_seeds", "pooled_rms_coarse", "pooled_rms_fine",
              "pooled_rms_difference", "pooled_relative_rms_difference", "rms_fine_seed_sd", "min_seed_relative_rms",
              "max_seed_relative_rms"), summary)


def plot_induced(output, raw):
    frames = SNAPSHOTS[1:]
    fig, axes = plt.subplots(4, len(frames), sharex=True, sharey="row", figsize=(12, 9))
    for col, frame in enumerate(frames):
        samples = raw[1, :, 3, frame]
        mean, sd = samples.mean(axis=0), samples.std(axis=0, ddof=1)
        for component, ax in enumerate(axes[:, col]):
            band = ax.fill_between(ZETA, mean[:, component] - sd[:, component], mean[:, component] + sd[:, component],
                                   color="#a8cfdf", alpha=0.55, label="Mean ± seed SD")
            traces = ax.plot(ZETA, samples[:, :, component].T, color="#477989", alpha=0.28, lw=0.7)
            line, = ax.plot(ZETA, mean[:, component], color="#123f59", lw=1.7, label="Ensemble mean")
            ax.axvline(0, color="0.52", lw=0.7, ls=":")
            ax.axhline(0, color="0.65", lw=0.6, zorder=0)
            ax.set_xlim(ZETA[0], ZETA[-1])
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            if col == 0:
                ax.set_ylabel(LABELS[component])
        axes[0, col].set_title(f"Physical t = {REPORTED[frame] + 0.004:.3f}")
    fig.suptitle("Grid B collision-induced EMT", x=0.1, ha="left", fontsize=17, weight="semibold")
    fig.text(0.1, 0.93, "Eight matched seeds at dT = 0.008, with both isolated controls subtracted", color="0.28")
    fig.supxlabel(r"Position relative to crossing $\zeta$ (simulation units)", y=0.062)
    fig.legend((line, band, traces[0]), ("Ensemble mean", "Mean ± seed SD", "Individual seeds"),
               loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=3, frameon=False)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.87, bottom=0.12, hspace=0.24, wspace=0.12)
    fig.savefig(output / "induced-profiles.png", dpi=170)
    plt.close(fig)


def plot_differences(output, time_map, delta):
    s = time_map[:, 4]
    seed_rms = rms(delta, axis=3)  # RMS over z within each seed, mode, time and component.
    colors = ("#295c9b", "#b56b00", "#247957", "#a23b76")
    fig, axes = plt.subplots(2, 2, sharex=True, figsize=(11, 7.4))
    for component, ax in enumerate(axes.flat):
        ax.axvline(0, color="0.55", lw=0.9, ls=":")
        for m, (mode, color) in enumerate(zip(MODES, colors)):
            samples = seed_rms[:, m, :, component]
            ax.fill_between(s, samples.min(axis=0), samples.max(axis=0), color=color, alpha=0.12, linewidth=0)
            ax.plot(s, samples.mean(axis=0), color=color, lw=1.55, label=mode.replace("-", " "))
        ax.set_title(LABELS[component], loc="left")
        ax.set_yscale("log")
        ax.set_xlim(s[0], s[-1])
        ax.grid(axis="y", which="major", color="0.9", lw=0.6)
        ax.set_axisbelow(True)
    fig.suptitle("Grid B paired timestep differences", x=0.095, ha="left", fontsize=17, weight="semibold")
    fig.text(0.095, 0.918, "Fine minus interpolated coarse, with the seed mean and range of the spatial RMS", color="0.28")
    fig.supylabel("RMS difference over longitudinal sites", x=0.019)
    fig.supxlabel(r"Time relative to nominal crossing $s$ (simulation units)", y=0.065)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=4, frameon=False)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.85, bottom=0.14, hspace=0.22, wspace=0.18)
    fig.savefig(output / "timestep-differences.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=WORKDIR / "runs" / "grid-b-ensemble-analysis",
                        help="fresh directory for the analysis outputs")
    output = parser.parse_args().output.resolve()
    if output.exists():
        parser.error(f"Output already exists: {output}. Choose a fresh --output directory.")

    queue = runpy.run_path(str(WORKDIR / "run-grid-b-ensemble.py"))
    print("Revalidating the archived Grid B runs", flush=True)
    completed = queue["completed_jobs"]()
    missing = set(queue["JOBS"]) - completed.keys()
    if missing:
        parser.error(f"Missing validated Grid B runs: {sorted(missing)}")
    output.mkdir(parents=True)
    print(f"Reading the {len(completed)} validated runs", flush=True)
    raw = load_campaign(queue, completed, output)
    time_map, coarse, fine = align_times(raw)
    delta = fine - coarse
    write_timestep_metrics(output, raw, time_map, coarse, fine)
    plot_induced(output, raw)
    plot_differences(output, time_map, delta)
    write_profiles(output, raw, time_map, coarse, fine, delta)
    print(f"Analysis written to {output}")


if __name__ == "__main__":
    main()
