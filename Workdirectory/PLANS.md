# Study plans and implementation handoff

5 November 2026

## Start here

The goal is a collision discretization and resolution-sensitivity demonstration using the existing SU(2) classical Yang–Mills code with minimal changes. The plans are ready to guide implementation and pilots. Their parameters are specified; numerical adequacy and scientific conclusions still require measurements.

| Document | Role and scope |
|---|---|
| [Collision plan with 2 to 1 spacing](plan-collision-spacing-2-to-1.md) | Active local series: 64 × 64 × 256, 96 × 96 × 384, 128 × 128 × 512; dT = 0.08; finest-grid dT = 0.04 check |
| [Collision plan with 4 to 1 spacing](plan-collision-spacing-4-to-1.md) | Corrected base series: 32 × 32 × 256 and 64 × 64 × 512 at dT = 0.08; Grid B ensemble follow-up is listed separately |
| [Grid B timestep and induced-EMT ensemble plan](plan-collision-spacing-4-to-1-grid-b-ensemble.md) | Completed 48-run Grid B study with eight seeds, dT = 0.08 and 0.008, and matched single-nucleus controls |
| [Roihu and seed 1 production handoff](ROIHU_REPRODUCTION.md#seed-1-production-handoff) | Historical solver comparison through the full-grid compatibility pilot |
| [Production lattice plan for seed 1](plan-production-seed1.md) | Completed `128 × 128 × 1024` triple at `dT=0.008` through `t=88` and [Grid B comparison](PRODUCTION_SEED1_HANDOFF.md) |
| [Single nucleus reference](plan-single-nucleus-reference.md) | Supporting diagnostics and historical study outline; not a required standalone campaign |

The ratio in each collision filename is a_perp:a_z. Both series keep the physical box (128, 128, 256). Do not add a 128 × 128 × 1024 run to the friend's series. Keep each series' time-step check distinct.

The corrected 4:1 checkpoint is summarized in [COLLISION_4_TO_1_HANDOFF.md](COLLISION_4_TO_1_HANDOFF.md). The full Grid B campaign and paired analysis are recorded in [GRID_B_ENSEMBLE_HANDOFF.md](GRID_B_ENSEMBLE_HANDOFF.md). The six pilot runs have a separate [seed-1 record](GRID_B_SEED1_HANDOFF.md).

## Shared decisions

- Preserve the established evolution equations, random-number generator and draw order, center formulas, boundaries, floating-point options and output columns. Use the continuum transverse-density convention: at fixed longitudinal profile and source strength, Var(rho_i) scales as 1/(a_x a_y), so the sampled density amplitude scales as 1/sqrt(a_x a_y). The first 4:1 batch used the opposite area factor and remains a pre-correction baseline; do not combine it with corrected runs.
- For both nuclei use Q_s = 0.125, Q_s σ_z = 0.5, σ_z = 4, m = 0.125 and Λ = 0.625 in simulation units. Q_s is the input source-strength convention; do not reinterpret 0.125 as Q_s a_perp or assign GeV/fm units without another specification.
- Retain i_L = int(0.85 N_z/2 − 1) and i_R = int(1.15 N_z/2 − 1), including truncation. Record physical positions, nominal crossing time and midpoint. Use raw and appropriately sampled aligned comparisons; a fixed-center rewrite is outside the current implementation scope.
- Keep the legacy generator. Identical seed and draw order match sources at the same transverse grid for collision/control and time-step comparisons. The same seed across different transverse meshes does not match physical realizations. A new Fourier-mode generator is not required or scheduled.
- Default to seed 1 for new pilots. Use another seed when explicitly matching documented historical data. Record the seed in every run. The new Grid B follow-up uses seeds 1–8, paired across time steps and run types as specified in its plan.
- Initially use eight CPU threads on the local i7-13700K. Run configurations sequentially. Benchmark the friend's own resources before choosing its fixed thread count.
- All full study runs, including the single-nucleus reference, use t_end = 88 from initialization. This gives Q_s t_end = 11 and nominally 1,100 steps at dT = 0.08. For collisions, Q_s times the elapsed time after nominal impact ranges from 8.5625 to 8.625 across these grids. The short coarse collision pilot still ends at t = 40, nominally 500 steps; its first 20–50 steps provide timing data.
- Collision fields are the first target. To isolate collision-induced T00, T0z and pressures, reuse matched left-only and right-only controls and subtract their tensors. This is supporting collision work, not a separate convergence campaign that must finish first.

## Current implementation state

The default build still selects `RightNucleusSimulation.cpp` in `SourceCode/SUNcSolver.cpp`. A build with `COLLISION_STUDY_BUILD` selects the existing collision procedure in `SourceCode/Simulation.cpp`; exactly one procedure include is active. An executable's name does not select its physics.

The collision-study build path in `SourceCode/LATTICE/Parameters.cpp` fixes `Qs0A`, `Qs0A_L`, and `Qs0A_R` to 0.125, `MEff` to 0.125, and `Lambda` to 0.625. Both initialization calls receive `Qs0A_L`; the `-QsA` option changes `Qs0A`, so the study launcher passes it explicitly and the configuration keeps the actual initializer amplitudes consistent. The standard build retains the existing defaults. Study executables require an explicit `COLLISION_STUDY_END_T` at compile time.

`SourceCode/COLLIDINGOBJECTS/CollidingObjects.cpp` now samples the transverse charge density with a `1/sqrt(a_x*a_y)` factor. The Poisson solve uses physical lattice momenta and consumes this array as a density. For a fixed longitudinal profile and nucleus multiplier, `(a_x*a_y)*Var(rho_i)` should equal `(Qs0A*T(z))^2` per color component. The completed 4:1 runs and analysis listed below predate this correction; retain them as historical baseline and rebuild before any corrected-source run.

A GSL check using seed 1, three color components per site, and raw transverse source draws before the longitudinal profile and smoothing measured `(a_x*a_y)*Var(rho_i)` = 0.0158164474 for Grid A (3,072 samples) and 0.0158337893 for Grid B (12,288 samples), against the target `Qs^2 = 0.015625` on both. Applying the old `sqrt(a_x*a_y)` factor to those same draws gives 4.04901054 and 0.253340629, respectively, against targets 4 and 0.25. Both executables were rebuilt with the corrected source, then the corrected pilot and two main runs completed as recorded below.

Existing runtime options include dimensions, spacings, `-dT`, both width inputs, cutoffs, seed, density multipliers and output directory. There is no endpoint option. The build script sets `EndSimulationT` separately for the pilot and main executables. The stopping loop is unchanged; a console marker reports actual final steps and time so the runner can record them. Floating-point accumulation can produce an extra final step, so compare common saved times. The 4:1 EMT profile writer now targets a physical cadence of 0.8 units (10 steps at dT=0.08; 100 at dT=0.008) and rejects time steps that cannot represent it within tolerance. Existing checkpoint outputs were produced before this cadence change.

CPU threading is already implemented and tested. `run-small.sh` honors `OMP_NUM_THREADS`; the solver uses OpenMP and threaded FFTW. A 100-step coarse-grid performance sweep took 200.9 seconds with one thread, 33.4 with eight and 30.4 with sixteen. Numerical comparisons passed. Results are in `runs/thread-benchmark-2c2gsgcw/`. These runs used the current short-test physics and dT = 0.008, not the study configuration.

`run-small.sh` still passes the original 16 × 16 × 256 smoke-test arguments. `benchmark-threads.py` expects source ending time 0.8. Use `build-collision-spacing-4-to-1.sh` to produce separate t=40 and t=88 collision executables, then `run-collision-spacing-4-to-1.sh A pilot`, `... A main`, and `... B main` to create fresh runs with manifests. The runner defaults to seed 1, eight threads, and dT=0.08; set `COLLISION_DT=0.008` for the fine-step case. Grid B main also accepts `collision`, `left-only`, and `right-only` run types; the eight-seed, two-time-step matrix is specified in [the new ensemble plan](plan-collision-spacing-4-to-1-grid-b-ensemble.md). Existing binaries and historical output must not be assumed to match the plans.

### Completed 4:1 initial batch

The Grid A t=40 pilot and both planned t=88 collision runs completed on 4 November 2026 with seed 1 and eight threads, before the transverse source normalization correction. Each manifest records the source revision and hashes, build recipe, actual endpoint and step count, resources, and complete solver command. Their results remain useful as the original batch baseline but must not be mixed with runs using the corrected density convention.

| Run | Steps / actual endpoint | Wall time | Peak resident memory |
|---|---:|---:|---:|
| Grid A pilot, t_end=40 | 501 / 40.0799999999995 | 26 s | 268,188 KiB |
| Grid A main, t_end=88 | 1,101 / 88.0799999999985 | 56 s | 268,184 KiB |
| Grid B main, t_end=88 | 1,101 / 88.0799999999985 | 552 s | 2,080,100 KiB |

The unchanged stopping loop took one additional step on each run. All main-run data files were finite; every Grid A profile has 256 rows and every Grid B profile has 512 rows. The main-run Gauss-law maximum/RMS residual peaks were 3.50845e-8 / 1.06815e-10 on Grid A and 4.54265e-11 / 5.60797e-13 on Grid B. These are diagnostics only, not an accuracy criterion.

Main outputs are under `runs/collision-spacing-4-to-1-A-main-t88-mOlDCj/` and `runs/collision-spacing-4-to-1-B-main-t88-ILFQYU/`; the pilot is under `runs/collision-spacing-4-to-1-A-pilot-t40-vbik8S/`. `analyze-collision-spacing-4-to-1.py` generated raw and aligned EMT profiles, crossing-plane histories, and Grid B minus Grid A data under `runs/collision-spacing-4-to-1-analysis/`. Recreate them with `python3 analyze-collision-spacing-4-to-1.py --grid-a runs/collision-spacing-4-to-1-A-main-t88-mOlDCj --grid-b runs/collision-spacing-4-to-1-B-main-t88-ILFQYU --output runs/collision-spacing-4-to-1-analysis`. The CSV compares equal reported profile times, aligns longitudinal coordinates, and records each grid's separate aligned time. It does not interpolate in time across the eight-unit profile interval.

This is a one-seed cross-grid comparison with different transverse meshes, so the differences combine sampling and resolution sensitivity and do not establish convergence. The optional Grid B dT=0.008 run was not part of this initial batch. Matched single-nucleus controls were not run, so these outputs do not isolate an induced EMT tensor.

### Corrected 4:1 rerun

The corrected density convention was run on 4 November 2026 with seed 1 and eight threads. Each fresh manifest records the corrected `CollidingObjects.cpp` hash and source convention.

| Run | Steps / actual endpoint | Wall time | Peak resident memory |
|---|---:|---:|---:|
| Grid A pilot, t_end=40 | 501 / 40.0799999999995 | 27 s | 268,000 KiB |
| Grid A main, t_end=88 | 1,101 / 88.0799999999985 | 58 s | 268,016 KiB |
| Grid B main, t_end=88 | 1,101 / 88.0799999999985 | 561 s | 2,080,268 KiB |

All runs completed with exit status 0 and finite output. The pilot has six 256-row profiles; each main run has twelve profiles, with 256 rows on Grid A and 512 rows on Grid B. Maximum Gauss-law residuals over time were 5.91042e-12 (RMS 1.00487e-13) for the pilot, 7.57964e-11 (RMS 1.25535e-12) for Grid A main, and 1.23941e-11 (RMS 1.41527e-13) for Grid B main.

Outputs are under `runs/collision-spacing-4-to-1-A-pilot-t40-aTVfkb/`, `runs/collision-spacing-4-to-1-A-main-t88-6gwTqr/`, and `runs/collision-spacing-4-to-1-B-main-t88-asxl1e/`.

The corrected pair was analyzed at 12 reported profile times in `runs/collision-spacing-4-to-1-corrected-analysis/`. The bulk site sums multiplied by cell volume give integrated energy 10.5014 (A) and 9.45867 (B) at t=0, a B/A ratio of 0.9007; at t=88 the ratio is 0.8863. This is far closer than the pre-correction initial ratio of 17.665. In aligned longitudinal profiles, the Grid B minus Grid A RMS difference relative to Grid A at t=24/40/88 was 10.9/11.6/14.2% for T00, 28.7/28.2/29.1% for T0Z, 6.4/17.9/9.4% for transverse pressure P_T, and 10.9/11.6/14.2% for longitudinal pressure P_L. This is a one-seed sensitivity comparison: different transverse meshes do not yield matched physical realizations, so it does not establish convergence. The later Grid B ensemble provides matched controls on Grid B.

## First implementation steps

1. Read the selected collision plan and [the source and output guide](GETTING_STARTED.md). Preserve existing local changes and historical data. Record the source revision and local configuration differences used for each executable.
2. The collision-study build configuration and thin launchers are prepared. They select `Simulation.cpp`, set the actual and recorded source amplitudes to 0.125, and build separate t=40 and t=88 executables. The only measurement addition is a console marker for the actual final step and time; the center, evolution, and output-sampling algorithms are unchanged.
3. Build or rebuild both executables from Ubuntu/WSL with `bash build-collision-spacing-4-to-1.sh` in `Workdirectory`. This compiles only `SourceCode/SUNcSolver.cpp` using the working recipe below.
4. The Grid A pilot was run with `OMP_NUM_THREADS=8 bash run-collision-spacing-4-to-1.sh A pilot`; its fresh output directory contains the manifest, log, resource report, and profiles. Unknown solver options are silently ignored, and `--help` is not a harmless help mode; verify the manifest and recorded values for future runs.
5. The pilot review found finite outputs, 256-row profiles, incoming propagation toward the center and post-impact central fields. It recorded 501 steps, actual t=40.0799999999995, 26 seconds, and 268,188 KiB peak RSS. The main-run Gauss-law and profile checks are summarized above. Small Gauss-law residuals alone do not establish physical accuracy.
6. The initial full-collision pair completed with `bash run-collision-spacing-4-to-1.sh A main` and `bash run-collision-spacing-4-to-1.sh B main`. The main executable uses `EndSimulationT = 88`. The optional `COLLISION_DT=0.008` Grid B extension remains separate from the initial batch.

After the source configuration in step 2 is prepared, a local 2:1 coarse pilot can use:

```bash
mkdir -p build runs
g++ -std=gnu++14 -O2 -fopenmp \
  SourceCode/SUNcSolver.cpp \
  -o build/CollSUNcSolver-pilot-t40 \
  -lgsl -lgslcblas -lfftw3_omp -lfftw3 -lm

export OMP_NUM_THREADS=8
export OMP_DYNAMIC=FALSE
run_dir=$(mktemp -d runs/collision-spacing-2-to-1-pilot-XXXXXX)
./build/CollSUNcSolver-pilot-t40 \
  -Nx 64 -Ny 64 -Nz 256 \
  -ax 2 -ay 2 -az 1 \
  -dT 0.08 -RlOverAzL 4 -RlOverAzR 4 \
  -QsA 0.125 -MEff 0.125 -Lambda 0.625 \
  -NucleusL 1 -NucleusR 1 -SEED 1 -o "$run_dir" \
  > "$run_dir/run.log" 2>&1
```

For the friend's 4:1 pilot, use its Grid A instead: dimensions 32, 32, 256 and spacings 4, 4, 1. The widths remain 4, 4. Neither command fixes the amplitude/endpoint wiring by itself; the source preparation and fresh build are essential.

## Outputs and completion of the first stage

Retain a manifest with the plan name, procedure, actual source parameters, dimensions/spacings, centers, seed, dT, endpoint, thread count, source/build identity and command. Keep each run in its own directory. `InfoID` contains only part of the required information.

The longitudinal EMT columns are z_index, T00, TXX, TYY, TZZ, T0X, T0Y, T0Z, followed by six alternative `mod` definitions. Use the first definitions consistently unless a separate comparison is intended. Bulk energy has a six-column initial row and four-column subsequent rows; parse these separately. Bulk T00 is a site sum, so normalize by site count for mean density or multiply by cell volume for integrated energy. The profile already averages over the transverse grid. Convert z_index to physical z.

The EMT profile writer uses `UOld` and `E`: after an evolution step, `UOld` is saved at the midpoint between link updates, so the evolved profile's physical measurement time is its reported time plus `dT/2`. The initial frame mixes initial links with an electric field at `dT/2` and is handled separately. The seed-1 comparison used 109 common physical times from `1.604` through `88.004`, linearly interpolating coarse profiles between 0.8-unit cadence brackets. Retain raw reported timestamps alongside any corrected physical-time estimate. Bulk diagnostics have their own staggering and must not automatically inherit the profile convention.

The first stage is complete when a correctly configured coarse collision has reached t = 40 and its manifest, profiles, diagnostic histories and measured costs have been reviewed. The full study additionally requires the main grid comparisons and its selected temporal check, plus matched controls when reporting the induced tensor.

The completed Grid B ensemble quantifies paired timestep differences and induced EMT uncertainty on the fixed Grid B mesh. Report one-seed comparisons between different transverse meshes as combined stochastic and spatial-resolution sensitivity. [Grid B results and analysis files](GRID_B_ENSEMBLE_HANDOFF.md)
