# Production lattice plan for seed 1

Run a matched collision, left-only and right-only triple through nominal `t=88` on the `128 × 128 × 1024` lattice. Use the existing seed-1 `64 × 64 × 512` triple at `dT=0.008` as the baseline. The optimized EMT writer passed the performance pilot below. The queue and cross-grid analysis are specified for implementation as one unattended job.

## Run configuration

| Setting | Production lattice | Grid B baseline |
| --- | --- | --- |
| Dimensions | `128 × 128 × 1024` | `64 × 64 × 512` |
| Spacings `(ax, ay, az)` | `(1, 1, 0.25)` | `(2, 2, 0.5)` |
| Box lengths | `128 × 128 × 256` | `128 × 128 × 256` |
| Longitudinal width in z cells, both nuclei | `16` | `8` |
| Longitudinal width in simulation units | `4` | `4` |
| Seed, timestep, nominal end time | `1`, `0.008`, `88` | `1`, `0.008`, `88` |
| Nominal crossing `(t, z)` | `(19.125, 127.625)` | `(19.25, 127.25)` |

Use the corrected SU(2) collision-study solver with `Qs=0.125`, `MEff=0.125` and `Lambda=0.625`, periodic x/y and open z. Run all three modes through `Simulation.cpp`. Set `-NucleusL/-NucleusR` to `(1,1)` for collision, `(1,0)` for left-only and `(0,1)` for right-only. Both source arrays are drawn in the same order in every mode, so the active nucleus in each control shares its realization with the collision on that mesh. [Grid B seed 1 handoff](GRID_B_SEED1_HANDOFF.md) · [Production pilot](ROIHU_REPRODUCTION.md#completed-full-grid-pilot)

The archived Grid B executable at `build/CollSUNcSolver-4to1-main-t88` matches the Grid B manifests. Its source tree at `2feb0d629954280be377836d9d1359ce0c2e276e` also matches the pre-optimization production pilot. Build a separate production executable from a pinned revision containing the optimized EMT writer, using `COLLISION_STUDY_BUILD` and `COLLISION_STUDY_END_T=88`. Record the revision, compiler, command, binary hash and `ldd` output, and use that executable for all three modes. [Build recipe](build-collision-spacing-4-to-1.sh)

## Performance pilot

The original production pilot spent about 483 of 710 seconds on the evolved profile in its 100-step block. The EMT writer now refreshes the z ghost planes once per profile, before the OpenMP loop. The gauge links are read-only during the measurement. The original writer issued four 1.5 MiB copies per site on this grid. [EMT writer](SourceCode/OBSERVABLES/EnergyMomentumTensor.cpp) · [Boundary copies](SourceCode/LATTICE/GaugeLinks.cpp) · [Pilot record](ROIHU_REPRODUCTION.md#completed-full-grid-pilot)

The optimized production-grid collision pilot used seed 1, eight threads and the archived pilot's toolchain and settings through `t=0.8`. Frame time was measured from the final Gauss sample to the completed EMT file.

| Measurement | Archived pilot | Optimized pilot |
| --- | ---: | ---: |
| Evolved EMT frame | 483.31 s | 5.76 s |
| Total wall time | 20m08.51s | 4m22.31s |
| Peak RSS | 15.79 GiB | 15.79 GiB |

Both production-grid EMT frames matched to roundoff, with a maximum absolute difference of `1.02e-19`. Bulk and Gauss-law histories matched at the written precision. All three Grid B modes also passed comparison with the archived coarse-step triple at `dT=0.08` through `t=24`. Across the 31 frames, the largest absolute change in the four collision-minus-controls profiles was `1.12e-19`. Use the optimized writer for the production build. The local validation record is in `runs/emt-writer-validation-mw766swy/`.

## Unattended execution

The optimized 100-step block with EMT output took 239.3 seconds. Extrapolating the collision pilot to 110 blocks per mode gives **about 7h20m per mode, or 22 hours serially**. Reserve **24–30 hours** of host uptime for this estimate. Run one solver at a time with `OMP_NUM_THREADS=8`, `OMP_DYNAMIC=FALSE`, `LC_ALL=C` and GSL `mt19937` seeded via `-SEED 1`. Require at least 19 GiB available inside WSL before each mode and 1 GiB free on the output volume.

1. Adapt the existing Grid B queue and validator for the production geometry. Give the proposed `run-production-seed1.py` a read-only `--dry-run` that checks the pinned build, dependencies, resources, analysis imports and the three exact Grid B inputs. Print the solver arguments and pending modes in order. [Grid B queue](run-grid-b-ensemble.py)
2. Have `--run` repeat the preflight under a campaign lock that the solver inherits. Run collision, left-only, then right-only. Create a fresh attempt directory under `runs/production-seed1-dt0008/` for each mode and capture the exact command, grid and source settings, revision, binary hash, start and finish times, exit status, actual endpoint, wall time, peak RSS and solver log. Write progress and the current attempt path to `campaign.log` as the queue advances.
3. Validate each run before advancing. On a failure, record the check and attempt path and stop. Rerunning the same command revalidates and skips completed modes. It starts an interrupted mode at `t=0` in a fresh directory because the normal solver path has checkpoint writes disabled. After the third validation, run the comparison automatically and record its outcome separately so an analysis retry can reuse all three runs. [Restart behavior](GETTING_STARTED.md)
4. Exercise the dry-run, duplicate-launch guard and resume behavior on a short test before the production launch. On launch day, keep Windows on AC and awake with Ubuntu running. Before leaving, confirm the background queue survives closing the launch terminal, reaches the first evolved profile at `t=0.8` and has memory headroom without active swapping.

Once the queue is implemented, use the path to `Workdirectory` inside the default Ubuntu WSL distribution.

```powershell
$productionPath = '/mnt/c/path/to/Paskahila/Workdirectory'
wsl.exe --cd $productionPath --exec python3 -B run-production-seed1.py --dry-run
if ($LASTEXITCODE -ne 0) { throw 'Production queue preflight failed.' }
Start-Process -FilePath wsl.exe -WindowStyle Hidden `
  -ArgumentList '--cd', ('"{0}"' -f $productionPath), '--exec', 'python3', '-u', '-B', 'run-production-seed1.py', '--run'
```

## Completion checks

- Require exit status `0`, `#COMPLETED ID=1`, eight threads, working FFTW threading, matching manifest and `InfoID1.txt` settings, and a finite endpoint consistent with `11,000` or `11,001` steps. Reject `#ERROR` markers. Record the actual end time. Floating-point accumulation can add one final step.
- Require all 111 EMT frames at reported times `0, 0.8, …, 88`, each with 1,024 ordered z rows and 14 finite columns. Check the initial longitudinal Gauss profile and the 551 bulk energy and Gauss-law samples at 20-step intervals. Record peak Gauss max and RMS and their times for each mode.

## Compare the two lattices

The comparison reads these three validated Grid B directories under `Workdirectory/runs/`.

| Mode | Grid B input |
| --- | --- |
| Collision | `collision-spacing-4-to-1-B-main-t88-7T9qSW` |
| Left-only | `collision-spacing-4-to-1-B-left-only-main-t88-e0Sqoy` |
| Right-only | `collision-spacing-4-to-1-B-right-only-main-t88-JjgmG6` |

1. On each native mesh, form `T_induced = T_collision − T_left-only − T_right-only` at each saved frame. Compare the full collision, the two controls and the induced profiles for `T00`, conventional `T0Z`, `P_T=(TXX+TYY)/2` and `P_L=TZZ`. Compare integrated bulk energy on its saved time axis after multiplying the bulk site sum by cell volume (`0.25` for production and `2` for Grid B). [Existing EMT reader and subtraction](analyze-grid-b-ensemble.py)
2. Overlay native profiles at their saved times and positions. For quantitative differences, align `ζ = z_index × az − z_crossing` and use overlap-weighted production averages in the fully covered Grid B bins. Compare evolved frames at `s = t_reported + 0.004 − t_crossing`. Use Grid B's evolved timestamps as targets and linearly interpolate the production frames within their common time range. This gives 109 matched times from Grid B's reported `t=1.6` to `88`. Keep the mixed-time `t=0` frame separate. [Time staggering and alignment](ROIHU_REPRODUCTION.md#compare-the-production-grid-fork-with-grid-b)
3. Save overlays for the initial state, the impact brackets and the matched times around `24`, `40` and `88`. Write per-frame absolute RMS, maximum absolute difference and relative RMS (`RMS(Δ)/RMS(Grid B)`) to CSV. Retain the denominator and flag near-zero reference profiles. Include the run audit and coordinate mapping in `runs/production-seed1-dt0008/analysis/`.

The same seed labels different transverse charge realizations on these two meshes. The comparison measures the combined effect of changing the mesh and the realization for this pair of events.
