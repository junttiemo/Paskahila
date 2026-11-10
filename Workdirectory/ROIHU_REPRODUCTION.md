# Roihu comparison and seed 1 handoff

Historical comparison record from 5 November 2026. The follow-up completed a matched production-grid triple and its comparison with Grid B using the optimized EMT writer. [Production seed 1 results](PRODUCTION_SEED1_HANDOFF.md)

A fresh build of `a9e81dc` served as the historical reference. The original Roihu output is unavailable. The seed-0 `timetest008Coll` and `timetest008Coll2` archives first appeared at `853d302` and provide a separate check at the historical transverse resolution.

## Measured agreement

The seed-0 collision pairs at `a9e81dc` and `bf15cbc` used the same GCC 15.2/OpenMP build and library environment in Ubuntu under WSL. The reduced grid was `16 × 16 × 512`, with one thread and the historical spacings and time step. Both runs produced 81 bulk samples and 17 EMT profiles through reported `t=12.8`, past nominal impact at `t≈9.625`. Bulk energy, bulk source, initial bulk stresses, and the ten EMT columns with unchanged definitions agreed exactly at the written precision. Both runs had peak Gauss-law max/RMS of `8.28e-10 / 1.16e-11`.

| Changed estimator | Relative RMS across the 17 profiles | Maximum absolute difference |
| --- | ---: | ---: |
| `T0Xmod` | 26.7% | `1.36e-8` |
| `T0Ymod` | 14.9% | `5.37e-8` |
| `T0Zmod` | 0.0136% | `8.29e-9` |

Each run took about 3m36s. The manifests, logs, raw data and per-frame comparisons are in `runs/roihu-pair-w925bujz/`.

The archival-grid pair used `128 × 128 × 512`, eight threads and a shortened endpoint of `t=0.8`. The bulk and Gauss-law files were identical, and the ten unchanged EMT columns differed by at most `1.11e-15` relative RMS across the two frames. The relative RMS differences for `T0Xmod`, `T0Ymod` and `T0Zmod` were 0.504%, 0.0330% and 0.000438%, respectively. The runs took 10m47s and 11m00s and both peaked at 7.91 GiB RSS. See `runs/roihu-pair-h3tmoi18/differences.csv`.

As a separate check, the rebuilt `a9e81dc` solver reproduced the earlier `timetest008Coll` archive with relative RMS differences of `1.34e-11` in the six bulk energy samples and `1.11e-11` in the two longitudinal energy-density profiles. See `runs/roihu-pair-h3tmoi18/archive-vs-reference.csv`.

## Historical configuration

The job generator at `a9e81dc:Workdirectory/Script.sh` specifies `128 × 128 × 1024`, spacings `(1, 1, 0.25)`, seeds 1–8 and separate collision, left-nucleus and right-nucleus executables. The source supplies SU(2), `dT=0.008`, longitudinal widths of 16 z cells (4 simulation units), source amplitude 0.125, infrared cutoff 0.125, ultraviolet cutoff 0.625 and an endpoint of 88. The script requests 40 CPUs per task. These are settings reconstructed from the commit.

The comparison used the fork's collision-study build, which selects the collision procedure and fixes the source amplitude at 0.125. Both collision initializers receive `Qs0A_L`. The runtime `-QsA` option sets `Qs0A`, a separate variable. [Build parameters](SourceCode/LATTICE/Parameters.cpp) · [Collision initialization](SourceCode/Simulation.cpp)

The two transverse charge-density factors, `sqrt(ax*ay)` and `1/sqrt(ax*ay)`, both equal 1 at the historical spacings. Commits `3ea8e64` and `8dac9e0` change the transported plaquettes used to compute `T0Xmod`, `T0Ymod` and `T0Zmod`. The comparisons record absolute and relative RMS differences separately for those three estimators. The relative RMS denominator is the reference RMS for the same component.

## Seed 1 production handoff

### Completed full-grid pilot

The pilot used the `128 × 128 × 1024` production lattice at spacings `(1, 1, 0.25)`, seed 1 and `dT=0.008`. Both solvers ran sequentially in Ubuntu under WSL 2 with GCC 15.2 and eight OpenMP threads. Each completed 100 steps through `t=0.8` and produced six bulk samples and two EMT profiles with 1024 longitudinal rows each.

| Solver | Wall time | Peak RSS | Initialization and initial profile | 100 steps and evolved profile |
|---|---:|---:|---:|---:|
| Reference `a9e81dc4467dd864c2e68803f166e2ddd01f27cf` | 20m08.33s | 15.79 GiB | 8m19s | 11m47s |
| Fork `04bd349b12a0661995a7696a00e8685d5f157e7a` | 20m08.51s | 15.79 GiB | 8m20s | 11m50s |

The phase timings and memory headroom come from three-second samples. Ubuntu had at least 4.40 GiB available in those samples, with zero pages swapped in or out. Each EMT profile took about eight minutes to compute.

The bulk energy, source term, initial bulk stresses and Gauss-law histories were identical. Both had peak Gauss-law max/RMS of `6.76e-14 / 4.95e-16`. The ten unchanged EMT columns agreed within `1.13e-15` relative RMS across the two profiles.

| Changed estimator | Relative RMS | Maximum absolute difference |
|---|---:|---:|
| `T0Xmod` | 63.1% | `6.99e-13` |
| `T0Ymod` | 7.04% | `6.65e-14` |
| `T0Zmod` | 0.000416% | `1.95e-10` |

The production-grid check covers initialization and early evolution before nominal impact at `t≈19.125`. The reduced-grid pair supplies the comparison through and after its impact.

The local, Git-ignored record is `runs/roihu-pair-sn3lhjkk/`. It contains `differences.csv`, `pilot-telemetry.tsv` and a `reference/` and `fork/` directory with `manifest.json`, `resource-usage.txt`, logs and raw outputs. The manifests preserve the exact commands, source revisions, compiler, binary hashes and linked libraries.

### Pre-optimization collision recipe through `t=88`

The original collision-only recipe pins the fork to `04bd349b12a0661995a7696a00e8685d5f157e7a`. Its `SourceCode` tree matches the completed Grid B run at `2feb0d629954280be377836d9d1359ce0c2e276e`. It uses the collision-study build with both nucleus multipliers set to 1, source amplitude 0.125, `σ_z/a_z=16`, `m=0.125` and `Λ=0.625`.

The pre-optimization estimate was roughly 22 hours for one solver, based on initialization plus 110 blocks of 100 steps and an EMT profile. The optimized timing is recorded with the new campaign. [EMT performance pilot](plan-production-seed1.md#performance-pilot)

From `Workdirectory` in Ubuntu, use Bash, Git, tar, g++ with OpenMP, GSL, FFTW and GNU time. This recipe builds the pinned source and records the commands, toolchain and outputs in a fresh run directory.

```bash
(
    set -euo pipefail
    revision=04bd349b12a0661995a7696a00e8685d5f157e7a
    mkdir -p runs
    run_dir=$(mktemp -d "$PWD/runs/roihu-fork-t88-XXXXXX")
    printf '%s\n' "$run_dir"
    exec 3>"$run_dir/commands.log"
    BASH_XTRACEFD=3
    set -x

    git -C .. archive "$revision" Workdirectory/SourceCode | tar -x --strip-components=1 -C "$run_dir"
    cd "$run_dir"
    printf 'source_revision=%s\n' "$revision" > MANIFEST.txt
    g++ --version >> MANIFEST.txt
    g++ -std=gnu++14 -O2 -fopenmp -DCOLLISION_STUDY_BUILD -DCOLLISION_STUDY_END_T=88 \
        SourceCode/SUNcSolver.cpp -o solver -lgsl -lgslcblas -lfftw3_omp -lfftw3 -lm > build.log 2>&1
    sha256sum solver >> MANIFEST.txt
    ldd ./solver >> MANIFEST.txt

    export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=1 LC_ALL=C
    /usr/bin/time -v -o resource-usage.txt ./solver \
        -Nx 128 -Ny 128 -Nz 1024 -ax 1 -ay 1 -az 0.25 -dT 0.008 \
        -RlOverAzL 16 -RlOverAzR 16 -QsA 0.125 -MEff 0.125 -Lambda 0.625 \
        -SEED 1 -NucleusL 1 -NucleusR 1 -o . > run.log 2>&1
)
```

Check `#COMPLETED ID=1`, the `#EVOLUTION_END_STEPS` marker, finite values and the Gauss-law history. Expect 111 EMT profiles of 1024 z rows each and 551 bulk samples through step 11000. Floating-point time accumulation can add a final evolution step. Record the actual endpoint and use the saved frames through nominal `t=88` for analysis.

### Compare the production-grid fork with Grid B

Reuse the completed seed-1 collision in `runs/collision-spacing-4-to-1-B-main-t88-7T9qSW/`. Its `MANIFEST.txt` and `InfoID1.txt` identify source revision `2feb0d629954280be377836d9d1359ce0c2e276e`, grid `64 × 64 × 512`, spacings `(2, 2, 0.5)`, `dT=0.008` and profiles through `t=88`. The physical box, source strength, cutoffs and physical widths match the production setup. [Grid B seed-1 record](GRID_B_SEED1_HANDOFF.md)

| Geometry | Grid B | Production grid |
|---|---:|---:|
| Nominal crossing time | 19.25 | 19.125 |
| Crossing position | 127.25 | 127.625 |
| Cell volume | 2 | 0.25 |

Compare the fork runs in simulation units. Convert the bulk site sum to integrated energy by multiplying by the cell volume. Plot `T00`, conventional `T0Z`, `P_T=(TXX+TYY)/2` and `P_L=TZZ` first at their saved coordinates and then on common coordinates relative to the crossing. Convert each z index to distance using its grid spacing and subtract that grid's crossing position. For evolved profiles, use time relative to impact as `t_reported + 0.004 − t_crossing`. Compare the initial frames separately. Use spatial averages over common coarse z bins and interpolate in time within the shared range. Inspect initialization, the incoming profiles, the closest frames around impact, and later profiles around `t=24, 40, 88`. The 0.8-unit output cadence limits resolution of differences right at impact.

Report absolute differences alongside relative RMS on the aligned grid. The same seed labels different transverse charge realizations on the two meshes, and both transverse and longitudinal spacings change. Treat the result as a comparison of these two events that combines mesh and realization effects.

The completed production campaign compares the collision-minus-controls observable using left-only and right-only runs on each lattice. [Cross-grid results](PRODUCTION_SEED1_HANDOFF.md#longitudinal-emt-profiles)
