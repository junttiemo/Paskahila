# Production lattice seed 1 results

The production triple and Grid B comparison finished at 03:04 Helsinki time on 9 November 2026, 22 h 12 min after launch. All three runs used the `128 × 128 × 1024` lattice at `dT=0.008` with eight OpenMP threads. The Grid B baseline used `64 × 64 × 512` at the same timestep and physical box size. [Production plan](plan-production-seed1.md) · [Grid B seed 1 record](GRID_B_SEED1_HANDOFF.md)

## Run audit

| Mode | Local attempt directory under `runs/production-seed1-dt0008/` | Wall time | Peak RSS | Peak Gauss max / RMS |
| --- | --- | ---: | ---: | ---: |
| Collision | `collision-rp664ymo` | 7 h 27 min 16 s | 15.79 GiB | `1.258e-13 / 2.313e-15` |
| Left-only | `left-only-dn7zxnhz` | 7 h 23 min 05 s | 15.79 GiB | `1.031e-13 / 7.568e-16` |
| Right-only | `right-only-cy51x8et` | 7 h 21 min 32 s | 15.79 GiB | `1.033e-13 / 7.147e-16` |

Each run exited with status 0 after 11,001 steps at `t≈88.008`. The queue validated 111 finite EMT frames of 1,024 rows and 14 columns, the initial longitudinal Gauss profile, and 551 bulk-energy and Gauss-law samples per mode. The two isolated controls' integrated energies stayed within 0.004% of their initial values. [Run manifests index](results/production-seed1/manifest-index.csv) · [Bulk energy histories](results/production-seed1/integrated-energy.csv)

The production executable came from a campaign-local Git snapshot of the Grid B source tree with the validated EMT writer applied. The source revision is `1109c22f6cace11ef7e6ac95c17fdc1a76315f69`. The build record contains the base revision, source tree, GCC 15.2 compile command, executable SHA-256, and linked libraries. The short queue test exercised duplicate-launch rejection, lock inheritance by the solver, and restarting an interrupted mode after skipping a validated collision. [Build record](results/production-seed1/build.json) · [Queue](run-production-seed1.py)

## Integrated energy

Each lattice's bulk site sum was multiplied by its cell volume, 0.25 for production and 2 for Grid B. The induced energy is the collision energy minus both isolated controls on the same lattice and saved bulk time.

| Run | Initial collision | Collision at t=88 | Sum of controls at t=88 | Induced at t=88 | Induced / sum of controls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Production seed 1 | 8.27657 | 8.68368 | 8.27684 | 0.40684 | 4.92% |
| Grid B seed 1 | 9.45875 | 9.70314 | 9.45999 | 0.24316 | 2.57% |

Production's induced energy was 67% higher than the Grid B seed 1 result. Across all eight fine-step Grid B triples, the induced energy ranged from 0.21253 to 0.42144, with a mean of 0.32643. The production value falls near the upper end of that range. The normalized production value of 4.92% is slightly above the observed Grid B range of 2.40–4.69%. The endpoint table uses the step-11,000 row of each archived `EnergyID<seed>.txt` and identifies each triple by grid and seed. [Energy plot](results/production-seed1/integrated-energy.png) · [Endpoint energies by seed](results/production-seed1/endpoint-energy-by-seed.csv) · [Grid B run index](results/grid-b-ensemble/manifest-index.csv)

## Longitudinal EMT profiles

The comparison forms the collision-minus-controls tensor on each native mesh. It compares the conventional `T00`, `T0Z`, `P_T=(TXX+TYY)/2`, and `P_L=TZZ`. Production profiles are averaged into the 511 fully covered Grid B bins after shifting `ζ=z−z_crossing`. Evolved profiles are matched at `s=t_reported+dT/2−t_crossing` by interpolating between the production frames bracketing each Grid B target. There are 109 evolved matches from Grid B reported `t=1.6` through `t=88`, plus a separate mixed-time initial comparison. [Comparison settings](results/production-seed1/comparison.json) · [Time map](results/production-seed1/time-mapping.csv) · [Spatial weights](results/production-seed1/z-mapping.csv)

At the final matched frame, both grids show two outgoing energy peaks and positive transverse pressure between them. The full collision `T00` differs by 12.4% RMS relative to Grid B seed 1. The absolute RMS difference in induced `T00` is 1.3% of the full Grid B collision `T00` RMS.

| Induced component at Grid B t=88 | Grid B seed 1 RMS | Absolute RMS difference | Difference / Grid B seed 1 RMS | Squared difference near outgoing fronts |
| --- | ---: | ---: | ---: | ---: |
| `T00` | `1.04e-7` | `1.07e-7` | 103% | 93.3% |
| `T0Z` | `8.50e-8` | `1.04e-7` | 122% | 95.0% |
| `P_T` | `2.76e-8` | `8.83e-9` | 32.1% | 14.3% |
| `P_L` | `7.79e-8` | `1.00e-7` | 129% | 95.8% |

The last column sums the squared difference in bins within six simulation units of either nominal outgoing front at `ζ=±68.754` and divides by the sum over all matched bins. For the ensemble comparison, subtract the pointwise mean of the eight fine-step Grid B profiles from the aligned production profile and take the spatial RMS. Dividing this by the RMS of the pointwise Grid B sample standard deviation gives 0.76, 0.75, 1.04 and 0.74 for `T00`, `T0Z`, `P_T` and `P_L`, respectively. Both calculations use the same 511 bins. [Aligned final profiles](results/production-seed1/aligned-t88.png) · [Per-frame RMS and maximum differences](results/production-seed1/profile-differences.csv)

The qualitative collision response agrees across the two runs. The final induced differences are comparable to the Grid B seed spread, which is much larger than the paired Grid B timestep differences at late times. Seed 1 generates different transverse charge fields on the two meshes, so this comparison measures mesh and realization effects together. [Grid B ensemble analysis](GRID_B_ENSEMBLE_HANDOFF.md#paired-results)

## Reproduction and archives

The plots and CSVs under [`results/production-seed1/`](results/production-seed1/) preserve the run audit, energy histories, alignment maps, and all 110 frame comparisons. Per-frame RMS ratios include the Grid B denominator and a flag for near-zero references. Use the absolute RMS for flagged rows.

The complete EMT outputs, Gauss-law histories, manifests, logs, and all snapshot overlays are in the local `runs/production-seed1-dt0008/` archive. The analysis accepts explicit production and Grid B directories in collision, left-only, right-only order. From `Workdirectory`, use Python with NumPy and Matplotlib and supply an existing empty output directory.

```bash
mkdir -p runs/production-seed1-recheck
python3 -B analyze-production-seed1.py \
  --output runs/production-seed1-recheck \
  --production runs/production-seed1-dt0008/collision-rp664ymo \
               runs/production-seed1-dt0008/left-only-dn7zxnhz \
               runs/production-seed1-dt0008/right-only-cy51x8et \
  --grid-b runs/collision-spacing-4-to-1-B-main-t88-7T9qSW \
           runs/collision-spacing-4-to-1-B-left-only-main-t88-e0Sqoy \
           runs/collision-spacing-4-to-1-B-right-only-main-t88-JjgmG6
```

To repeat the solver campaign, restore the three archived Grid B inputs listed in the plan and follow its pinned build and WSL launch recipe. Running the queue again with the validated production archive reuses those attempts. [Production launch and resume recipe](plan-production-seed1.md#unattended-execution)
