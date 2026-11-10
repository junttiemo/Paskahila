# Grid B timestep and collision-induced EMT ensemble plan

5 November 2026

All 48 runs and the paired analysis are complete as of 7 November 2026. Results and output paths are in [GRID_B_ENSEMBLE_HANDOFF.md](GRID_B_ENSEMBLE_HANDOFF.md). This follow-up studies the corrected 4:1 Grid B at two time steps and uses matched single-nucleus controls to define the collision-induced energy-momentum tensor. The [seed-1 pilot](GRID_B_SEED1_HANDOFF.md), [earlier two-grid checkpoint](COLLISION_4_TO_1_HANDOFF.md), and [shared study plan](PLANS.md) record the implementation history.

## Questions and estimands

1. How much do the longitudinal EMT profiles at `dT=0.08` and `dT=0.008` differ on the same Grid B source realization?
2. What is the operational collision-induced profile, defined componentwise as
   `T_induced = T_collision - T_left-only - T_right-only`?
3. How much do these paired quantities vary across eight independent legacy-generator seeds?

The subtraction defines the requested comparison; because the field evolution is nonlinear, it is not assumed to be a unique decomposition of the full interacting field. Report it as a collision-minus-isolated-controls observable.

## Fixed configuration

| Setting | Value |
|---|---:|
| Gauge group | SU(2) |
| Grid | `64 x 64 x 512` |
| Spacings `(a_x,a_y,a_z)` | `(2,2,0.5)` |
| Physical box `(L_x,L_y,L_z)` | `(128,128,256)` |
| Source strength for each active nucleus | `Q_s=0.125` simulation units |
| Width inputs `sigma_z/a_z` | `8` for both nuclei (`sigma_z=4`) |
| Infrared and ultraviolet regulators | `m=0.125`, `Lambda=0.625` |
| End time | `t_end=88` from initialization |
| Time steps | `dT=0.08` and `dT=0.008` |
| Seeds | `1,2,3,4,5,6,7,8` |
| Threads | Eight, fixed for the campaign |
| EMT profile cadence | `0.8` simulation-time units |

For the 4:1 implementation keep the accepted center formulas and boundary treatment. On Grid B the center indices are `(216,293)`, the nominal physical centers are `(108,146.5)`, the aligned midpoint is `z_c=127.25`, and the nominal crossing time is `t_coll=19.25`. Keep raw `t,z` coordinates and also record `s=t-t_coll` and `zeta=z-z_c` for aligned analysis.

The corrected source convention is cell-centered transverse density: the Gaussian amplitude scales as `1/sqrt(a_x*a_y)`, so `(a_x*a_y)*Var(rho_i)` equals `(NucleusMultiplier*Q_s*T(z))^2` before smoothing. Confirm the compiled source and manifest use this convention before launching. The current tracked corrected-code results in the handoff used the inverse-area normalization; do not mix the earlier pre-correction batch into this ensemble.

## Run matrix and matching

Run all combinations below for every seed:

| Time step | Collision | Left-only | Right-only | Runs across 8 seeds |
|---|---:|---:|---:|---:|
| `0.08` | 8 | 8 | 8 | 24 |
| `0.008` | 8 | 8 | 8 | 24 |
| **Total** | **16** | **16** | **16** | **48** |

For each `(seed,dT)` triple, use exactly the same seed and grid in all three run types. Set density multipliers to `(1,1)` for collision, `(1,0)` for left-only, and `(0,1)` for right-only. The existing initializer draws both transverse source arrays in the same order even when one multiplier is zero, so each active nucleus retains the corresponding realization from the collision run. The source draws are independent of `dT`, so reusing a seed also pairs the source realization across the two time steps. Preserve this code path and verify the multipliers in each manifest.

The existing launcher supports these modes on Grid B main runs:

```bash
COLLISION_SEED=1 COLLISION_DT=0.08 bash run-collision-spacing-4-to-1.sh B main collision
COLLISION_SEED=1 COLLISION_DT=0.08 bash run-collision-spacing-4-to-1.sh B main left-only
COLLISION_SEED=1 COLLISION_DT=0.08 bash run-collision-spacing-4-to-1.sh B main right-only
```

Repeat for seeds 1 through 8, then repeat the same triples with `COLLISION_DT=0.008`. Run sequentially at eight threads to stay within the measured memory footprint. Keep every run in its own fresh directory and retain its command, revision and source hashes, seed, multiplier pair, actual endpoint/step count, elapsed time and peak RSS in `MANIFEST.txt`.

### Unattended continuation for seeds 2–8

The ensemble queue requires the archived seed-1 executable, its build record, all six seed-1 run directories, and a checkout with the same committed `SourceCode` tree as `2feb0d6`. It runs pending coarse triples followed by fine triples and checks each run's manifest, profile cadence, finite output, and diagnostic histories before advancing. It appends progress and Gauss-law peaks to `runs/grid-b-ensemble.log`. Restarting revalidates completed runs and creates a fresh directory for each pending run. [Ensemble queue](run-grid-b-ensemble.py)

To use the queue from Windows PowerShell, set the path to your `Workdirectory` as it appears inside WSL. Preview the queue before launching the background process.

```powershell
$gridBPath = '/mnt/c/path/to/Paskahila/Workdirectory'
wsl.exe --cd $gridBPath --exec python3 -B run-grid-b-ensemble.py --dry-run
Start-Process -FilePath wsl.exe -WindowStyle Hidden `
  -ArgumentList '--cd', $gridBPath, '--exec', 'python3', '-u', 'run-grid-b-ensemble.py', '--run'
```

The commands use the default WSL distribution. Keep Windows awake and Ubuntu running for about 46 hours. The queue runs independently of the PowerShell and Codex sessions. Run the same start command to resume after an interruption. The campaign log lists each validated run directory and the final completion marker.

## Output schedule and time alignment

Save the EMT profile every `0.8` simulation-time units: every 10 steps at `dT=0.08`, and every 100 steps at `dT=0.008`. The initial profile at `t=0` is retained. Both resolutions should therefore provide the same nominal profile times from `t=0` through `t=88`, 111 profiles per run or 5,328 profiles in the campaign. The existing 20-step bulk-energy and Gauss-law outputs keep their step-based schedule; they are diagnostics, not the EMT analysis series.

The evolution loop can take one extra final step because of floating-point accumulation. Record the actual marker from each log; the primary profile comparison uses the common saved frames through nominal `t=88`, not an assumption about the final loop time. After an evolution step, the EMT profile's physical measurement time is its reported time plus `dT/2`; the initial frame mixes times and is treated separately. The seed-1 comparison used 109 common physical times from `1.604` through `88.004`, linearly interpolating coarse profiles between 0.8-unit cadence brackets. Retain raw reported timestamps alongside any corrected physical-time estimate.

The seed-1 block with the planned 111-profile cadence took about 38 minutes per run at `dT=0.08` and 93 minutes per run at `dT=0.008`, with eight threads and about 1.98 GiB peak RSS. The observed averages gave a planning estimate of 46 hours for the 42-run continuation. The queue completed that block in 43 h 24 min. The [ensemble handoff](GRID_B_ENSEMBLE_HANDOFF.md#run-audit) records the actual runtimes and memory by run type.

## Analysis

For each seed and time step, calculate `T_induced` componentwise for `T00`, `T0Z`, `P_T=(TXX+TYY)/2`, and `P_L=TZZ`. Convert the longitudinal index to physical `z`, and produce both raw and center-aligned profiles. Compare the two time steps within each run type and for `T_induced` at matched physical times. Keep the full EMT tensor columns in archived raw files.

Across seeds, report the ensemble mean and standard deviation or standard error for each paired timestep difference and for `T_induced`. Use seed-paired statistics: form each seed's control subtraction and fine-minus-coarse difference before aggregating. Show individual seed results as well as the aggregate so a small mean does not hide seed variation. Do not present the eight seeds as eight independent samples for each of the three terms in a control subtraction; each seed contributes one matched triple.

Review the initial profiles and Gauss-law residuals for each run, check all output values are finite, and verify profiles occur at the requested common physical times. Report measured runtime and memory by time step and run type. The existing corrected one-seed comparison is a checkpoint, not an ensemble result; this plan does not by itself establish spatial convergence or a universal physical prediction.

## Execution gates and deliverables

1. The seed-1 build passed the corrected source and cadence gates. The queue verified and reused the same executable for seeds 2–8.
2. All 16 matched triples passed the multiplier, cadence, finite-output, and completion checks.
3. The [ensemble handoff](GRID_B_ENSEMBLE_HANDOFF.md) links the 48-run manifest index, profile-time audit, paired raw and aligned CSVs, ensemble plots, and uncertainty summaries.
