#!/usr/bin/env python3
"""Validate and resume the Grid B campaign using the seed-1 executable."""

import argparse
import hashlib
from itertools import product
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


WORKDIR = Path(__file__).resolve().parent
RUNS = WORKDIR / "runs"
LAUNCHER = "run-collision-spacing-4-to-1.sh"
BUILD = "build/CollSUNcSolver-4to1-main-t88"
PILOT_REVISION = "2feb0d629954280be377836d9d1359ce0c2e276e"
PILOT_BINARY_SHA256 = "09818024fe227a1d60195555421f67c1d037da28d3884f381c439892c42f8594"
BUILD_SOURCES = (
    "SourceCode/SUNcSolver.cpp",
    "SourceCode/LATTICE/Parameters.cpp",
    "SourceCode/Simulation.cpp",
    "SourceCode/COLLIDINGOBJECTS/CollidingObjects.cpp",
)
MULTIPLIERS = {"collision": "1 1", "left-only": "1 0", "right-only": "0 1"}
TIME_STEPS = {"0.08": (1100, 10), "0.008": (11000, 100)}
JOBS = tuple(product(TIME_STEPS, map(str, range(1, 9)), MULTIPLIERS))
RUN_PATTERN = "collision-spacing-4-to-1-B-*-t88-*"
LOG_FORMAT = "%(asctime)s %(message)s"


def read_fields(path):
    return dict(line.split("=", 1) for line in path.read_text(encoding="utf-8").splitlines() if "=" in line)


def require_fields(fields, expected):
    for key, value in expected.items():
        if fields.get(key) != value:
            raise ValueError(f"{key} is {fields.get(key)!r}, expected {value!r}")


def verify_build():
    record = WORKDIR / f"{BUILD}.build-info.txt"
    require_fields(read_fields(record), {
        "role": "main-t88",
        "endpoint": "88.0",
        "selected_procedure": "Simulation.cpp",
        "source_revision": PILOT_REVISION,
        "source_density_convention": "cell-centered transverse density",
    })
    hashes = {path: digest for digest, path in re.findall(
        r"^([0-9a-f]{64})  (\S+)$", record.read_text(encoding="utf-8"), re.MULTILINE
    )}
    if hashes.get(BUILD) != PILOT_BINARY_SHA256:
        raise ValueError("The build record does not identify the seed-1 executable")
    for path in (*BUILD_SOURCES, BUILD):
        if hashlib.sha256((WORKDIR / path).read_bytes()).hexdigest() != hashes.get(path):
            raise ValueError(f"{path} differs from its seed-1 build record hash")
    head_tree, pilot_tree = subprocess.check_output([
        "git", "rev-parse", "HEAD:Workdirectory/SourceCode", f"{PILOT_REVISION}:Workdirectory/SourceCode",
    ], cwd=WORKDIR, text=True).splitlines()
    if head_tree != pilot_tree:
        raise ValueError("The committed SourceCode tree differs from the seed-1 revision")
    if b"\r\n" in (WORKDIR / LAUNCHER).read_bytes():
        raise ValueError(f"{LAUNCHER} needs LF line endings for Bash")


def read_numeric_rows(path, width, initial_width=None):
    rows = []
    with path.open(encoding="utf-8") as stream:
        for row_number, line in enumerate(stream, 1):
            columns = initial_width if row_number == 1 and initial_width is not None else width
            values = tuple(map(float, line.split()))
            if len(values) != columns or not all(map(math.isfinite, values)):
                raise ValueError(f"{path.name}:{row_number} must have {columns} finite values")
            rows.append(values)
    return rows


def check_z_rows(path, width, nz=512):
    rows = read_numeric_rows(path, width)
    if len(rows) != nz or any(row[0] != index for index, row in enumerate(rows)):
        raise ValueError(f"{path.name} must contain z indices 0 through {nz - 1} in order")
    return rows


def validate_diagnostics(run_dir, seed, dt, nominal_steps, profile_steps, nz):
    profile_names = [f"EMTDynT{frame * 8 / 10:g}ID{seed}.txt"
                     for frame in range(nominal_steps // profile_steps + 1)]
    if {path.name for path in run_dir.glob("EMTDynT*")} != set(profile_names):
        raise ValueError(f"expected {len(profile_names)} EMT profiles at 0.8-unit intervals")
    initial = check_z_rows(run_dir / profile_names[0], 14, nz)
    for name in profile_names[1:]:
        check_z_rows(run_dir / name, 14, nz)
    check_z_rows(run_dir / f"GVDynT0ID{seed}.txt", 2, nz)

    gauss = read_numeric_rows(run_dir / f"GausCheck{seed}.txt", 3)
    energy = read_numeric_rows(run_dir / f"EnergyID{seed}.txt", 4, initial_width=6)
    if len(gauss) != nominal_steps // 20 + 1 or len(energy) != len(gauss):
        raise ValueError("bulk energy or Gauss-law history has the wrong length")
    for index, (gauss_row, energy_row) in enumerate(zip(gauss, energy)):
        time = 20 * index * float(dt)
        energy_time = energy_row[1] if index else energy_row[0]
        if energy_row[0] != 20 * index or not (
            math.isclose(gauss_row[0], time, abs_tol=1e-7)
            and math.isclose(energy_time, time, abs_tol=1e-7)
        ):
            raise ValueError(f"missing diagnostic sample at t={time:g}")
        if not 0 <= gauss_row[2] <= gauss_row[1]:
            raise ValueError(f"invalid Gauss-law residuals at t={time:g}")
    return initial, gauss


def validate_run(run_dir, job, manifest):
    dt, seed, mode = job
    nominal_steps, profile_steps = TIME_STEPS[dt]
    require_fields(manifest, {
        "run_status": "completed", "exit_status": "0", "completion_marker": "yes",
        "procedure": "Simulation.cpp", "stage": "main", "grid": "Grid B", "model": "SU(2)",
        "executable_hash": f"{PILOT_BINARY_SHA256}  {BUILD}",
        "source_revision": PILOT_REVISION,
        "source_density_convention": "cell-centered transverse density",
        "Nx_Ny_Nz": "64 64 512", "ax_ay_az": "2 2 0.5", "sigma_z_over_az_left_right": "8 8",
        "center_indices_left_right": "216 293", "endpoint_from_initialization": "88",
        "aligned_center_z_c": "127.25", "nominal_collision_time": "19.25",
        "actual_initializer_amplitudes_left_right": "0.125 0.125 (both calls use Qs0A_L)",
        "m": "0.125", "Lambda": "0.625",
        "seed": seed, "dT": dt, "run_type": mode,
        "charge_density_multipliers_left_right": MULTIPLIERS[mode],
        "threads": "8", "OMP_DYNAMIC": "FALSE", "emt_profile_cadence_simulation_time": "0.8",
        "emt_profile_cadence_steps": str(profile_steps), "nominal_steps": str(nominal_steps),
    })
    steps = int(manifest.get("actual_steps", "0"))
    endpoint = manifest.get("actual_endpoint", "nan")
    if steps not in (nominal_steps, nominal_steps + 1) or not math.isclose(
        float(endpoint), steps * float(dt), abs_tol=1e-7
    ):
        raise ValueError(f"unexpected endpoint {endpoint} at step {steps}")
    seconds = int(manifest.get("elapsed_seconds", "-1"))
    rss = int(manifest.get("peak_rss_kb", "0"))
    if seconds < 0 or rss <= 0:
        raise ValueError("missing wall time or peak RSS")

    log = set((run_dir / "run.log").read_text(encoding="utf-8").splitlines())
    markers = {
        "#NUMBER OF THREADS 8", "#FFTW THREAD STATUS 1", f"#COMPLETED ID={seed}",
        f"#EMT_PROFILE_CADENCE_T=0.8 STEPS={profile_steps}", f"#EVOLUTION_END_STEPS={steps} TIME={endpoint}",
    }
    if not markers.issubset(log) or any(line.startswith("#ERROR") for line in log):
        raise ValueError("run.log lacks expected markers or reports an error")
    require_fields(read_fields(run_dir / f"InfoID{seed}.txt"), {
        "#Nx,Ny,Nz": " 64 64 512", "#ax,ay,az": " 2 2 0.5", "#dT": f" {dt}",
    })

    initial, gauss = validate_diagnostics(run_dir, seed, dt, nominal_steps, profile_steps, 512)

    logging.info(
        "Validated %s | seed %s dT=%s %s | 111 profiles | %ss, %s KiB RSS | "
        "initial peak |T00| %.6g | peak Gauss max/RMS %.6g / %.6g",
        run_dir.name, seed, dt, mode, seconds, rss, max(abs(row[1]) for row in initial),
        max(row[1] for row in gauss), max(row[2] for row in gauss),
    )


def completed_jobs():
    completed = {}
    for run_dir in sorted(RUNS.glob(RUN_PATTERN)):
        path = run_dir / "MANIFEST.txt"
        if not path.is_file():
            continue
        try:
            manifest = read_fields(path)
            job = tuple(manifest.get(key) for key in ("dT", "seed", "run_type"))
            if job not in JOBS or manifest.get("run_status") != "completed":
                continue
            validate_run(run_dir, job, manifest)
        except (OSError, ValueError) as error:
            logging.warning("Ignoring %s: %s", run_dir.name, error)
            continue
        if job in completed:
            raise ValueError(f"Multiple validated runs for {job}: {completed[job].name}, {run_dir.name}")
        completed[job] = run_dir
    return completed


def launch_job(job, lock_fd):
    dt, seed, mode = job
    env = dict(os.environ, COLLISION_SEED=seed, COLLISION_DT=dt, OMP_NUM_THREADS="8", OMP_DYNAMIC="FALSE", LC_ALL="C")
    env.pop("GSL_RNG_TYPE", None)
    env.pop("GSL_RNG_SEED", None)
    run_dir = None
    with subprocess.Popen(
        ["bash", LAUNCHER, "B", "main", mode], cwd=WORKDIR, env=env, pass_fds=(lock_fd,),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    ) as process:
        for line in process.stdout:
            logging.info("%s", line.rstrip())
            match = re.fullmatch(r"Writing (runs/collision-spacing-4-to-1-B-[^/]+) \(8 threads\)", line.strip())
            if match:
                run_dir = WORKDIR / match[1]
        if process.wait() != 0 or run_dir is None:
            raise ValueError(f"Launcher failed for seed {seed}, dT={dt}, {mode}. See the queue log and run directory.")
    validate_run(run_dir, job, read_fields(run_dir / "MANIFEST.txt"))


def run_queue(lock_fd=None):
    verify_build()
    completed = completed_jobs()
    if any(job[1] == "1" and job not in completed for job in JOBS):
        raise ValueError("All six archived seed-1 runs must pass validation before continuing")
    pending = [job for job in JOBS if job not in completed]
    logging.info("%s/48 runs validated, %s remaining | seed-1 executable SHA-256 %s",
                 len(completed), len(pending), PILOT_BINARY_SHA256)
    if os.name == "posix":
        if not os.access(WORKDIR / BUILD, os.X_OK) or not os.access("/usr/bin/time", os.X_OK):
            raise ValueError("The solver and /usr/bin/time must be executable in Ubuntu")
        libraries = subprocess.check_output(["ldd", BUILD], cwd=WORKDIR, text=True)
        if "not found" in libraries:
            raise ValueError(f"The seed-1 executable has missing shared libraries:\n{libraries}")
    if shutil.disk_usage(RUNS).free < len(pending) * 20 * 1024**2:
        raise ValueError("Insufficient free space for the remaining run outputs")
    for index, job in enumerate(pending, 1):
        dt, seed, mode = job
        logging.info("%s %s/%s | seed %s dT=%s %s",
                     "Starting" if lock_fd is not None else "Queued", index, len(pending), seed, dt, mode)
        if lock_fd is not None:
            verify_build()
            launch_job(job, lock_fd)
    if lock_fd is not None:
        logging.info("Grid B campaign complete. All 48 runs validated.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--dry-run", action="store_true", help="validate existing runs and list the queue (default)")
    modes.add_argument("--run", action="store_true", help="launch or resume the queue in Ubuntu under WSL")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, stream=sys.stdout)

    if not args.run:
        run_queue()
        return
    if os.name != "posix":
        parser.error("Use Ubuntu under WSL to run the solver")
    import fcntl

    with (RUNS / "grid-b-ensemble.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("A Grid B queue or one of its solver processes is already running") from error
        handler = logging.FileHandler(RUNS / "grid-b-ensemble.log", encoding="utf-8")
        handler.setFormatter(logging.Formatter(LOG_FORMAT))
        logging.getLogger().addHandler(handler)
        run_queue(lock.fileno())


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        logging.error("Grid B queue stopped: %s", error)
        sys.exit(1)
