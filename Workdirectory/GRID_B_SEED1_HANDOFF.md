# Grid B seed 1 handoff: timestep and control pilot

Pilot record from 5 November 2026. The completed Grid B timestep and collision-induced EMT ensemble has a separate [handoff](GRID_B_ENSEMBLE_HANDOFF.md).

## Status

The seed-1 block contains collision, left-only, and right-only runs at both `dT=0.08` and `dT=0.008`. All six runs used the corrected 4:1 collision-study source, seed 1, eight threads, and a 0.8-unit EMT profile cadence. They completed successfully, produced finite numeric outputs, and have matching source and build hashes. These six runs are included in the 48-run ensemble.

The raw output directories and sequence log are under the local, Git-ignored `Workdirectory/runs/` directory. The tracked record below preserves the run identities and key checks; the raw profiles are not part of this PR.

## Run record

| `dT` | Run type | Output directory | Steps / actual endpoint | Runtime | Peak RSS | Peak Gauss max / RMS |
|---:|---|---|---:|---:|---:|---:|
| 0.08 | Collision | `collision-spacing-4-to-1-B-main-t88-KtaLD9` | 1,101 / 88.0799999999985 | 38m23s | 2,080,324 KiB | 1.2394e-11 / 1.4153e-13 |
| 0.08 | Left-only | `collision-spacing-4-to-1-B-left-only-main-t88-xtgPRf` | 1,101 / 88.0799999999985 | 38m36s | 2,079,892 KiB | 1.2222e-11 / 1.3939e-13 |
| 0.08 | Right-only | `collision-spacing-4-to-1-B-right-only-main-t88-qK8bl6` | 1,101 / 88.0799999999985 | 38m25s | 2,080,256 KiB | 3.3243e-12 / 4.0754e-14 |
| 0.008 | Collision | `collision-spacing-4-to-1-B-main-t88-7T9qSW` | 11,001 / 88.0079999999947 | 1h32m31s | 2,080,032 KiB | 6.0975e-13 / 5.9100e-15 |
| 0.008 | Left-only | `collision-spacing-4-to-1-B-left-only-main-t88-e0Sqoy` | 11,001 / 88.0079999999947 | 1h33m10s | 2,080,284 KiB | 2.4544e-13 / 2.5750e-15 |
| 0.008 | Right-only | `collision-spacing-4-to-1-B-right-only-main-t88-JjgmG6` | 11,001 / 88.0079999999947 | 1h32m37s | 2,080,184 KiB | 5.3748e-13 / 5.2869e-15 |

Each run wrote 111 EMT profiles (512 rows × 14 columns) at the planned cadence, from the initial profile through nominal `t=88`. The 114 numeric output files in each run parsed as finite values. Logs had no error markers. The Gauss-law values are constraint diagnostics, not by themselves a physical accuracy test.

The fine-step sequence directory is `runs/grid-b-seed1-dt0008-triple-Vee48s/`. Its `SEQUENCE.txt` records `sequence_status=completed`, `exit_status=0`, and the order collision → left-only → right-only; `sequence.log` preserves the run output. The three coarse runs completed in the same order before the overnight sequence.

## Reading the measurements

The code's EMT profile time is staggered relative to the filename label. After an evolution step, `UOld` is saved at the midpoint between link updates; the EMT writer uses `UOld` and `E`. The initializer constructs those fields from links separated by one time step. Therefore an evolved saved profile's measurement time is its reported time plus `dT/2`. The initial profile at reported `t=0` mixes the initial link time and electric field at `dT/2`, so it is handled separately.

For the matched timestep comparison, 109 common physical times from `1.604` through `88.004` were compared. The coarse profiles were linearly interpolated in time between their 0.8-unit cadence brackets; the z grids are identical, so no spatial interpolation was needed. The earliest fine evolved frame at `0.804` has no earlier evolved coarse bracket and was omitted. The comparison retained the initial frame as a separate case.

Across the 109 matched times and all 512 z positions, RMS(fine − coarse) / RMS(coarse) for the operational induced tensor was:

| Component | Relative RMS difference |
|---|---:|
| `T00` | 0.14% |
| `T0Z` | 0.20% |
| `P_T=(TXX+TYY)/2` | 0.07% |
| `P_L=TZZ` | 0.22% |

Near the collision, at physical time `19.204`, the largest component snapshot difference was 1.13% for `T00`; `T0Z`, `P_T`, and `P_L` differed by 0.53%, 0.36%, and 0.27%. At `24.004`, `40.004`, and `88.004`, all four component differences were below 0.25%. This is strong timestep stability evidence for this one source realization, not a formal convergence result or ensemble uncertainty estimate.

The operational induced tensor is `T_collision − T_left-only − T_right-only`. At reported times 24, 40, and 88 for `dT=0.08`, the induced RMS values were:

| Component | RMS at 24 / 40 / 88 | Approximate share of collision-profile RMS |
|---|---:|---:|
| `T00` | 7.25e-8 / 1.16e-7 / 1.04e-7 | 0.9% / 1.4% / 1.3% |
| `T0Z` | 5.21e-8 / 7.15e-8 / 8.50e-8 | 0.8% / 1.0% / 1.2% |
| `P_T` | 7.56e-8 / 3.88e-8 / 2.75e-8 | Isolated-control `P_T` is near zero; subtraction leaves a small absolute signal |
| `P_L` | 1.01e-7 / 6.19e-8 / 7.78e-8 | 1.2% / 0.8% / 1.0% |

The subtraction is the study's operational collision-minus-isolated-controls observable. Since field evolution is nonlinear, it is not a unique decomposition of the interacting field. The induced profiles are small compared with the full collision profiles, so preserve absolute scales and inspect per-seed behavior in the remaining ensemble.

## Reproduction and continuation

For a reproduction of the Grid B pilot, build from the `2feb0d6` source tree with the collision-study flags. From `Workdirectory` in Ubuntu/WSL, run a coarse matched triple with the resulting executable.

```bash
for mode in collision left-only right-only; do
  COLLISION_SEED=1 COLLISION_DT=0.08 OMP_NUM_THREADS=8 \
    bash run-collision-spacing-4-to-1.sh B main "$mode"
done
```

The fine matched triple runner is [run-grid-b-seed1-dt0008-triple.sh](run-grid-b-seed1-dt0008-triple.sh). It checks the build record and source/binary hashes, creates a fresh sequence directory, runs the three modes in order, verifies each manifest and cadence, and stops if a check fails. The seed 1 fine sequence has already completed; do not rerun it just to continue the campaign.

The [ensemble queue](run-grid-b-ensemble.py) completed seeds 2–8 with the seed-1 build. The [ensemble handoff](GRID_B_ENSEMBLE_HANDOFF.md) records the measured runtimes, analysis methods, and results.

The detailed manifests, profile files, Gauss-law histories, and sequence logs are in the local `runs/` tree. See [PLANS.md](PLANS.md) for the broader study context and [the ensemble plan](plan-collision-spacing-4-to-1-grid-b-ensemble.md) for the run matrix and aggregation rules.
