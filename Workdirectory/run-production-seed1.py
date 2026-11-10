#!/usr/bin/env python3
"""Build, validate, and resume the seed 1 production lattice campaign."""

import argparse
from dataclasses import dataclass
from datetime import datetime
import hashlib
import io
import json
import logging
import math
import os
from pathlib import Path
import re
import runpy
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time


WORKDIR = Path(__file__).resolve().parent
RUNS = WORKDIR / "runs"
CAMPAIGN = RUNS / "production-seed1-dt0008"
GRID_B = runpy.run_path(str(WORKDIR / "run-grid-b-ensemble.py"))
BASELINE = {
    "collision": RUNS / "collision-spacing-4-to-1-B-main-t88-7T9qSW",
    "left-only": RUNS / "collision-spacing-4-to-1-B-left-only-main-t88-e0Sqoy",
    "right-only": RUNS / "collision-spacing-4-to-1-B-right-only-main-t88-JjgmG6",
}
BASE_REVISION = GRID_B["PILOT_REVISION"]
EMT_SOURCE = Path("SourceCode/OBSERVABLES/EnergyMomentumTensor.cpp")
EMT_BLOB = "578e217182b48215579a737f527529848a4bd145"
ENVIRONMENT = {
    "OMP_NUM_THREADS": "8", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "8",
    "GSL_RNG_TYPE": "mt19937", "GSL_RNG_SEED": "1", "LC_ALL": "C",
}
DT = 0.008
PROFILE_STEPS = 100
GIB = 1024**3
ANALYZER = "analyze-production-seed1.py"
LOG_FORMAT = "%(asctime)s %(message)s"
ERRORS = (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError)


@dataclass(frozen=True)
class Lattice:
    dimensions: tuple
    spacings: tuple
    width: int
    end_time: float

    @property
    def centers(self):
        nz = self.dimensions[2]
        return int(0.85 * nz / 2 - 1), int(1.15 * nz / 2 - 1)

    @property
    def crossing(self):
        left, right = (i * self.spacings[2] for i in self.centers)
        return (right - left) / 2, (left + right) / 2

    @property
    def nominal_steps(self):
        return round(self.end_time / DT)

    def solver_args(self, mode):
        left, right = GRID_B["MULTIPLIERS"][mode].split()
        return [
            "-Nx", str(self.dimensions[0]), "-Ny", str(self.dimensions[1]), "-Nz", str(self.dimensions[2]),
            "-ax", f"{self.spacings[0]:g}", "-ay", f"{self.spacings[1]:g}", "-az", f"{self.spacings[2]:g}",
            "-dT", f"{DT:g}", "-RlOverAzL", str(self.width), "-RlOverAzR", str(self.width),
            "-QsA", "0.125", "-MEff", "0.125", "-Lambda", "0.625",
            "-SEED", "1", "-NucleusL", left, "-NucleusR", right, "-o", ".",
        ]


PRODUCTION = Lattice((128, 128, 1024), (1, 1, 0.25), 16, 88)
COARSE = Lattice((64, 64, 512), (2, 2, 0.5), 8, 88)


def read_json(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def write_json(path, data):
    pending = path.with_suffix(".tmp")
    pending.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    pending.replace(path)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_output(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def compiler_command(lattice):
    return [
        "g++", "-std=gnu++14", "-O2", "-fopenmp", "-DCOLLISION_STUDY_BUILD",
        f"-DCOLLISION_STUDY_END_T={lattice.end_time:g}", "SourceCode/SUNcSolver.cpp",
        "-o", "../solver", "-lgsl", "-lgslcblas", "-lfftw3_omp", "-lfftw3", "-lm",
    ]


def linked_libraries(binary):
    libraries = subprocess.check_output(["ldd", str(binary)], text=True)
    if "not found" in libraries:
        raise ValueError(f"Missing libraries for {binary}\n{libraries}")
    return libraries.splitlines()


def verify_source_revision(source, revision):
    baseline = git_output(WORKDIR.parent, "ls-tree", "-r", f"{BASE_REVISION}:Workdirectory", "--", "SourceCode")
    old_blob = git_output(WORKDIR.parent, "rev-parse", f"{BASE_REVISION}:Workdirectory/{EMT_SOURCE.as_posix()}")
    expected = baseline.replace(f"{old_blob}\t{EMT_SOURCE.as_posix()}", f"{EMT_BLOB}\t{EMT_SOURCE.as_posix()}")
    if git_output(source, "ls-tree", "-r", revision, "--", "SourceCode") != expected:
        raise ValueError("The pinned source differs from the Grid B revision with the validated EMT writer")


def verify_build(campaign, lattice):
    build_dir = campaign / "build"
    source = build_dir / "source"
    binary = build_dir / "solver"
    build = read_json(build_dir / "build.json")
    revision = build["source_revision"]
    GRID_B["require_fields"](build, {
        "base_revision": BASE_REVISION, "compile_command": compiler_command(lattice),
        "source_revision": git_output(source, "rev-parse", "HEAD"),
        "source_tree": git_output(source, "rev-parse", f"{revision}:SourceCode"),
        "binary_sha256": sha256(binary),
    })
    verify_source_revision(source, revision)
    if not os.access(binary, os.X_OK) or not os.access("/usr/bin/time", os.X_OK):
        raise ValueError("The solver and GNU time must be executable")
    linked_libraries(binary)
    return build


def prepare_build(campaign, lattice):
    if (campaign / "build").exists():
        verify_build(campaign, lattice)
        logging.info("Verified the existing production build")
        return
    optimized = (WORKDIR / EMT_SOURCE).read_bytes().replace(b"\r\n", b"\n")
    blob = subprocess.check_output(["git", "hash-object", "--stdin"], input=optimized).decode().strip()
    if blob != EMT_BLOB:
        raise ValueError(f"{EMT_SOURCE} differs from the validated optimized source")

    staging = Path(tempfile.mkdtemp(prefix="build-attempt-", dir=campaign))
    source = staging / "source"
    source.mkdir()
    logging.info("Preparing pinned source and executable in %s", staging)
    archive = subprocess.check_output(["git", "archive", f"{BASE_REVISION}:Workdirectory", "SourceCode"], cwd=WORKDIR.parent)
    with tarfile.open(fileobj=io.BytesIO(archive)) as bundle:
        bundle.extractall(source, filter="data")
    (source / EMT_SOURCE).write_bytes(optimized)

    # Keep the source revision in the campaign's own repository.
    git_output(source, "init", "--quiet", "--initial-branch=production-seed1")
    git_output(source, "-c", "core.autocrlf=false", "-c", "core.filemode=false", "add", "SourceCode")
    git_output(source, "-c", "user.name=Production build", "-c", "user.email=production@localhost",
               "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false", "commit", "--quiet",
               "-m", f"Pin seed 1 solver from {BASE_REVISION} with validated EMT writer")
    revision = git_output(source, "rev-parse", "HEAD")
    verify_source_revision(source, revision)

    command = compiler_command(lattice)
    env = dict(os.environ, LC_ALL="C")
    with (staging / "build.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = staging / "solver"
    build = {
        "base_revision": BASE_REVISION, "source_revision": revision,
        "source_tree": git_output(source, "rev-parse", f"{revision}:SourceCode"),
        "compiler": subprocess.check_output(["g++", "--version"], text=True, env=env).splitlines()[0],
        "compile_command": command, "compile_working_directory": "build/source",
        "binary_sha256": sha256(binary), "linked_libraries": linked_libraries(binary),
    }
    write_json(staging / "build.json", build)
    staging.rename(campaign / "build")
    logging.info("Built %s | revision %s | SHA-256 %s", binary.name, revision, build["binary_sha256"])


def check_resources(campaign, pending):
    free = shutil.disk_usage(campaign if campaign.exists() else RUNS).free
    if free < GIB:
        raise ValueError(f"Output volume has {free / GIB:.2f} GiB free; need at least 1 GiB")
    if pending:
        info = Path("/proc/meminfo").read_text(encoding="utf-8")
        available = int(re.search(r"^MemAvailable:\s+(\d+) kB", info, re.MULTILINE)[1]) * 1024
        if available < 19 * GIB:
            raise ValueError(f"Ubuntu has {available / GIB:.2f} GiB available; need at least 19 GiB")
        logging.info("Resources: %.2f GiB available in Ubuntu, %.2f GiB free on output volume", available / GIB, free / GIB)


def validate_info(attempt, lattice):
    GRID_B["require_fields"](GRID_B["read_fields"](attempt / "InfoID1.txt"), {
        "#Nx,Ny,Nz": " " + " ".join(map(str, lattice.dimensions)),
        "#ax,ay,az": " " + " ".join(f"{a:g}" for a in lattice.spacings),
        "#dT": f" {DT:g}", "QsR": " 0.5", "QsAperp": " 0.125",
        "MEff aka Infrared Cutoff": " 1 (times QsAperp)",
        "Lambda aka UV Cutoff": " 5 (times QsAperp)",
        "SigmaLOverAz L": f" {lattice.width}", "SigmaLOverAz R": f" {lattice.width}",
    })


def run_settings(campaign, lattice, build, mode):
    return {
        "mode": mode, "seed": 1, "procedure": "Simulation.cpp", "model": "SU(2)",
        "dimensions": list(lattice.dimensions), "spacings": list(lattice.spacings),
        "box_lengths": [n * a for n, a in zip(lattice.dimensions, lattice.spacings)],
        "source_width_in_cells": [lattice.width, lattice.width], "source_amplitudes": [0.125, 0.125],
        "source_multipliers": list(map(int, GRID_B["MULTIPLIERS"][mode].split())),
        "MEff": 0.125, "Lambda": 0.625, "boundaries": ["periodic", "periodic", "open"],
        "center_indices": list(lattice.centers), "nominal_crossing_t_z": list(lattice.crossing),
        "dT": DT, "nominal_end_time": lattice.end_time, "environment": ENVIRONMENT,
        "source_revision": build["source_revision"], "binary_sha256": build["binary_sha256"],
        "command": [str(campaign / "build" / "solver"), *lattice.solver_args(mode)],
    }


def validate_attempt(attempt, lattice, expected, manifest):
    GRID_B["require_fields"](manifest, {**expected, "exit_status": 0})
    steps, endpoint = manifest["actual_steps"], manifest["actual_endpoint"]
    if steps not in (lattice.nominal_steps, lattice.nominal_steps + 1) or not math.isclose(
        float(endpoint), steps * DT, rel_tol=0, abs_tol=1e-7
    ):
        raise ValueError(f"unexpected endpoint {endpoint} at step {steps}")
    if not math.isfinite(manifest["wall_seconds"]) or manifest["wall_seconds"] < 0 or manifest["peak_rss_kib"] <= 0:
        raise ValueError("missing wall time or peak RSS")

    markers = {
        "#NUMBER OF THREADS 8", "#FFTW THREAD STATUS 1", "#GAUGE GROUP IS SU(2)", "#COMPLETED ID=1",
        "GSL_RNG_TYPE=mt19937", "GSL_RNG_SEED=1", "#SETTING UP LATTICE WITH OPEN BOUNDARY CONDITIONS",
        f"#EMT_PROFILE_CADENCE_T=0.8 STEPS={PROFILE_STEPS}", f"#EVOLUTION_END_STEPS={steps} TIME={endpoint}",
    }
    log = set((attempt / "run.log").read_text(encoding="utf-8").splitlines())
    if markers - log:
        raise ValueError(f"missing run.log markers {sorted(markers - log)}")
    if any(line.startswith("#ERROR") for line in log):
        raise ValueError("run.log contains #ERROR")
    validate_info(attempt, lattice)
    _, gauss = GRID_B["validate_diagnostics"](attempt, "1", DT, lattice.nominal_steps, PROFILE_STEPS, lattice.dimensions[2])
    max_time, max_value, _ = max(gauss, key=lambda row: row[1])
    rms_time, _, rms_value = max(gauss, key=lambda row: row[2])
    return {"max": max_value, "max_time": max_time, "rms": rms_value, "rms_time": rms_time}


def completed_modes(campaign, lattice, build):
    completed = {}
    for mode in BASELINE:
        for path in sorted(campaign.glob(f"{mode}-*/manifest.json")):
            try:
                manifest = read_json(path)
                if manifest.get("status") != "validated":
                    continue
                peaks = validate_attempt(path.parent, lattice, run_settings(campaign, lattice, build, mode), manifest)
                GRID_B["require_fields"](manifest, {"gauss_peaks": peaks})
            except ERRORS as error:
                raise ValueError(f"Cannot revalidate {path.parent}: {error}") from error
            if mode in completed:
                raise ValueError(f"Multiple validated {mode} attempts: {completed[mode]} and {path.parent}")
            completed[mode] = path.parent
            logging.info("Revalidated %s; skipping %s", path.parent, mode)
    return completed


def preflight(campaign, lattice, analysis_python):
    build = verify_build(campaign, lattice)
    logging.info("Verified production revision %s | SHA-256 %s", build["source_revision"], build["binary_sha256"])
    subprocess.run([analysis_python, "-B", ANALYZER, "--check-imports"], cwd=WORKDIR, check=True)
    for mode, path in BASELINE.items():
        try:
            GRID_B["validate_run"](path, (f"{DT:g}", "1", mode), GRID_B["read_fields"](path / "MANIFEST.txt"))
            validate_info(path, COARSE)
        except ERRORS as error:
            raise ValueError(f"Grid B baseline {path}: {error}") from error
    completed = completed_modes(campaign, lattice, build)
    pending = [mode for mode in BASELINE if mode not in completed]
    check_resources(campaign, pending)
    logging.info("Pending modes in order: %s", ", ".join(pending) or "none")
    for mode in pending:
        settings = run_settings(campaign, lattice, build, mode)
        logging.info("%s | %s", mode, shlex.join(settings["command"]))
    return build, completed, pending


def now():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def launch_mode(campaign, lattice, build, mode, lock_fd):
    settings = run_settings(campaign, lattice, build, mode)
    if sha256(Path(settings["command"][0])) != build["binary_sha256"]:
        raise ValueError("Production binary changed after preflight")
    attempt = Path(tempfile.mkdtemp(prefix=f"{mode}-", dir=campaign))
    path = attempt / "manifest.json"
    manifest = dict(settings, started=now(), status="running")
    logging.info("Starting %s | %s", mode, attempt)
    write_json(path, manifest)
    start = time.monotonic()
    try:
        with (attempt / "run.log").open("wb") as log:
            result = subprocess.run(
                ["/usr/bin/time", "-v", "-o", "resource-usage.txt", *settings["command"]],
                cwd=attempt, env={**os.environ, **ENVIRONMENT}, stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT, pass_fds=(lock_fd,),
            )
        manifest.update(exit_status=result.returncode, finished=now(), wall_seconds=time.monotonic() - start)
        resources = (attempt / "resource-usage.txt").read_text(encoding="utf-8")
        rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", resources)
        if rss:
            manifest["peak_rss_kib"] = int(rss[1])
        end = re.search(r"^#EVOLUTION_END_STEPS=(\d+) TIME=(\S+)$",
                        (attempt / "run.log").read_text(encoding="utf-8"), re.MULTILINE)
        if end:
            manifest.update(actual_steps=int(end[1]), actual_endpoint=end[2])
        peaks = validate_attempt(attempt, lattice, settings, manifest)
    except ERRORS as error:
        manifest.update(status="failed", finished=now(), failure=str(error))
        write_json(path, manifest)
        raise ValueError(f"{mode} attempt {attempt}: {error}") from error
    manifest.update(status="validated", gauss_peaks=peaks)
    write_json(path, manifest)
    logging.info("Validated %s | t=%s, %s steps | %.1fs, %.2f GiB peak RSS | Gauss max/RMS %.6g / %.6g",
                 attempt, manifest["actual_endpoint"], manifest["actual_steps"], manifest["wall_seconds"],
                 manifest["peak_rss_kib"] / 1024**2, peaks["max"], peaks["rms"])
    return attempt


def analyze(campaign, completed, analysis_python):
    output = campaign / "analysis"
    inputs = {
        "production": [completed[mode].relative_to(WORKDIR).as_posix() for mode in BASELINE],
        "grid_b": [path.relative_to(WORKDIR).as_posix() for path in BASELINE.values()],
    }
    if output.exists():
        GRID_B["require_fields"](read_json(output / "analysis-run.json"), {"status": "completed", **inputs})
        logging.info("Comparison already completed | %s", output)
        return

    staging = Path(tempfile.mkdtemp(prefix="analysis-attempt-", dir=campaign))
    command = [
        analysis_python, "-B", ANALYZER, "--output", staging.relative_to(WORKDIR).as_posix(),
        "--production", *inputs["production"], "--grid-b", *inputs["grid_b"],
    ]
    status = dict(inputs, status="running", command=command, started=now())
    record = campaign / "analysis-status.json"
    write_json(record, dict(status, attempt=str(staging)))
    logging.info("Starting production vs. Grid B comparison | %s", staging)
    try:
        with (staging / "analysis.log").open("wb") as log:
            result = subprocess.run(command, cwd=WORKDIR, stdout=log, stderr=subprocess.STDOUT, check=True)
        status.update(status="completed", finished=now(), exit_status=result.returncode)
        write_json(staging / "analysis-run.json", status)
        staging.rename(output)
    except ERRORS as error:
        status.update(status="failed", finished=now(), failure=str(error))
        if isinstance(error, subprocess.CalledProcessError):
            status["exit_status"] = error.returncode
        write_json(record, dict(status, attempt=str(staging)))
        raise ValueError(f"Comparison failed; see {staging / 'analysis.log'}") from error
    write_json(record, dict(status, attempt=str(output)))
    logging.info("Comparison completed | %s", output)


def locked_campaign(campaign, lattice, analysis_python, prepare=False):
    import fcntl

    campaign.mkdir(parents=True, exist_ok=True)
    with (campaign / "campaign.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("A production queue or one of its solvers is already running") from error
        handler = logging.FileHandler(campaign / "campaign.log", encoding="utf-8")
        handler.setFormatter(logging.Formatter(LOG_FORMAT))
        logger = logging.getLogger()
        logger.addHandler(handler)
        try:
            if prepare:
                prepare_build(campaign, lattice)
            else:
                logging.info("Queue started | PID %s", os.getpid())
                build, completed, pending = preflight(campaign, lattice, analysis_python)
                for mode in pending:
                    check_resources(campaign, pending=True)
                    completed[mode] = launch_mode(campaign, lattice, build, mode, lock.fileno())
                logging.info("All three modes validated")
                analyze(campaign, completed, analysis_python)
        except ERRORS as error:
            logging.error("Stopped: %s", error)
            raise
        finally:
            logger.removeHandler(handler)
            handler.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--prepare", action="store_true", help="build the pinned production solver")
    actions.add_argument("--dry-run", action="store_true", help="check prerequisites and list pending modes")
    actions.add_argument("--run", action="store_true", help="run or resume the queue and comparison")
    parser.add_argument("--analysis-python", default=sys.executable, help="Python with NumPy and Matplotlib")
    args = parser.parse_args()
    if not sys.platform.startswith("linux"):
        parser.error("Run this queue in Ubuntu under WSL")
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, stream=sys.stdout)
    if args.dry_run:
        preflight(CAMPAIGN, PRODUCTION, args.analysis_python)
    else:
        locked_campaign(CAMPAIGN, PRODUCTION, args.analysis_python, args.prepare)


if __name__ == "__main__":
    try:
        main()
    except ERRORS as error:
        sys.exit(f"Production queue stopped: {error}")
