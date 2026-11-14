# LOVECLIM experiment: coast or herds with a modelled paleoclimate

Run 14 November 2026 with `code/run_loveclim.sh` (400 long + 400 short chronology draws, paired across strategies and climates).

## What changed

The original model takes habitat from hand-drawn biome boxes, one Saharo-Arabian humidity index and one monsoon index. Here, per-cell and per-1-kyr net primary production, desert fraction and temperature come from the LOVECLIM transient run behind Timmermann & Friedrich 2016 (*Nature*, doi:10.1038/nature19365), regridded from 1° to 20′ (`code/climate_ibs.py`). The climate now drives both uses of "biome" (`code/model_ibs.py`):

1. **Big-game habitat:** `0.02 + 0.88 · q_npp · q_dry · q_cold · (1 − forest_pen · forest)`.
2. **Re-tooling classes:** hyper-arid, semi-arid, savanna or closed forest, set per cell and per slice from the climate. A population pays `biome_pen` years when it moves into a new class.

Coastal fresh water in dry phases follows the local desert fraction. There are six new sampled parameters (`npp_ref`, `des_exp`, `t_min`, `npp_forest`, `forest_pen`, `cl_shift`), drawn from a separate RNG. Sea level, coastlines, rivers, crossings and the solver are unchanged. The hand-drawn (`box`) rows reproduce `results/monte_carlo_runs.csv` to a relative difference of 4e-5. The gap comes from the sea-level table extended below 35 ka (erectus patch), which a negative time jitter can reach; the evidence-fit numbers are identical.

**Earlier variant.** `habitat_only/` holds the first LOVECLIM run (13 November). It used LOVECLIM habitat but kept the hand-drawn re-tooling classes (`run(..., biome_mode='static')`). That mode reproduces its stored runs to 1e-16. Its headline numbers differ (see the last table), and the earlier page text quoted them.

## Results (full LOVECLIM, long chronology unless stated)

| | Hand-drawn | LOVECLIM |
|---|---|---|
| Coastal reaches Sahul before megafauna | 75% | **82.5%** |
| Coastal win rate by coastal-productivity quartile (low → high) | 36, 79, 90, 95% | 71, 82, 85, 92% |
| Coastal win rate by exit date (earliest → latest quartile) | 72, 67, 80, 81% | 57, 84, 90, 99% |
| Megafauna reaches Sahul | 99% | 90% |
| Megafauna Horn → Sahul transit, median | 13.3 kyr | 16.9 kyr |
| Interior Arabia occupied by 70 ka: coastal / megafauna / mixed | 59 / 95 / 96% | 11 / 26 / 33% |
| Short chronology: coastal before megafauna | 21.7% | 13.5% |

- **Shore richness matters less.** Even when the coast is poorest, coastal wins 71% of runs, against 36% on the hand-drawn climate.
- **Exit timing and the Arabian desert decide it.** LOVECLIM's Arabian belt is about 4–15% desert at 82–84 ka, 50–66% at 76–79 ka and 78–85% at 66–75 ka. Only early exits give herds a chance.
- **The Nefud constraint is hard to meet in every strategy.** Under LOVECLIM, interior Arabia is occupied by 70 ka in 11–33% of runs, against 59–96% on the hand-drawn climate.

## Review tests on LOVECLIM (`review_ibs.json`)

Evidence filter: Sulawesi by 67.8 ka, Sumatra by 63 ka and, unless stated, the Nefud occupied by 70 ka.

| Mixed cost `c_mix` | P(mixed first at Sahul) | Posterior coastal / megafauna / mixed | Same, without the Nefud test |
|---|---|---|---|
| 0 | 0.99 | 0.16 / 0.36 / 0.48 | 0.33 / 0.26 / 0.41 |
| 0.1 | 0.82 | 0.17 / 0.38 / 0.45 | 0.34 / 0.27 / 0.39 |
| 0.2 | 0.19 | 0.18 / 0.41 / 0.41 | 0.35 / 0.29 / 0.36 |

On the hand-drawn climate the matching posteriors are 0.23 / 0.33 / 0.44, and 0.31 / 0.30 / 0.39 without the Nefud test.

- **Mixed wins only when flexibility is free.** At a 20% speed cost it reaches Sahul first in 19% of runs, and its fit ties megafauna.
- **The Nefud test still decides coastal vs megafauna.** With it, coastal fits worst (0.16–0.18). Without it, coastal fits better than megafauna (0.33–0.35 against 0.26–0.29) and ties mixed at 20% cost.
- **Few runs are accepted.** With the Nefud test, 248 of 1,200 strategy runs pass at c_mix 0, against 572 on the hand-drawn climate.

## Habitat-only variant vs full LOVECLIM

| | Habitat only (13 Nov) | Full LOVECLIM |
|---|---|---|
| Coastal before megafauna, long | 91.7% | 82.5% |
| Coastal win, poorest-shore quartile | 82% | 71% |
| Short chronology, coastal before megafauna | 32.5% | 13.5% |
| Posterior coastal / megafauna / mixed | 0.16 / 0.35 / 0.48 | 0.16 / 0.36 / 0.48 |

## Limits

- **One coarse climate simulation, not an ensemble**, so climate uncertainty is not propagated beyond a ±2 kyr shift.
- **The productivity-to-habitat mapping is assumed.**
- **Some small Lesser Sunda islands are barren in LOVECLIM.** They are refilled from the nearest resolved land (`climate_ibs.py`).
- **The Timmermann & Friedrich benchmark (`summary_ibs.json`, `benchmark`) shares the same climate family**, so it is not independent confirmation.
- **Raw climate files are not in the repo** (CC BY-NC, 130 MB). See `data/external/IBS_HDS_README.md` to fetch them.
