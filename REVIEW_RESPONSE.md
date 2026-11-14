# Response to the external review (14 Nov 2026)

All changes are behind flags. Defaults reproduce the published results exactly: rerunning c_mix = 0 and wait = 0 matches `results/monte_carlo_runs.csv` on all 800 runs (worst relative difference 3.6e-14), and `exp.json` `post_strategy`/`abc` are identical. Runner: `code/run_review.sh`; summary: `results/review_summary.json`.

## 1. Mixed strategy wins by construction (fixed: cost added)

Mixed speed is now (1 − c_mix) × the pointwise maximum of the two pure surfaces (`--c-mix`).

| Setting | P(mixed reaches Sahul first) | ABC posterior, coastal / terrestrial / mixed | Mixed Sahul arrival, median ka |
|---|---|---|---|
| c_mix 0, wait 0 kyr | 0.79 | 0.23 / 0.33 / 0.43 | 59.0 |
| c_mix 0.05, wait 0 kyr | 0.79 | 0.24 / 0.34 / 0.42 | 59.2 |
| c_mix 0.1, wait 0 kyr | 0.74 | 0.24 / 0.34 / 0.41 | 59.5 |
| c_mix 0.2, wait 0 kyr | 0.40 | 0.25 / 0.36 / 0.39 | 58.0 |
| c_mix 0, wait 2 kyr | 0.79 | 0.23 / 0.33 / 0.44 | 59.0 |
| c_mix 0, wait 5 kyr | 0.79 | 0.23 / 0.33 / 0.44 | 59.0 |
| c_mix 0, wait 10 kyr | 0.79 | 0.23 / 0.33 / 0.44 | 59.0 |
| c_mix 0.1, wait 5 kyr | 0.74 | 0.24 / 0.34 / 0.41 | 59.5 |

At c_mix = 0, mixed arrives first in 79% of runs by construction. A 20% efficiency cost cuts that to 40%. The ABC posterior barely separates the strategies at any cost (mixed 0.39–0.44, terrestrial 0.33–0.36, coastal 0.23–0.25). **"Mixed wins" is not supported as a strong claim.** It holds at zero cost by construction, and the evidence ratio is weak at every cost.

## 2. Waiting at a coast for lower sea level (fixed)

`--wait W` lets a population that can survive at a coastal node depart at the lowest-sea-level slice within W kyr. It has **no material effect**: P(mixed first) moves 0.789 → 0.790 and the posteriors move by at most 0.001 for W = 2, 5, 10.

## 3. Nefud ≥ 70 ka constraint vs the 70–73 ka humid window (tested)

The ABC keeps draws with Nefud arrival ≥ 70 ka, and the humid-window end `h5a_b` ~ U(70, 73) gates that arrival, so they are not independent. Dropping the constraint (`exp.py TAG --no-nefud`) changes the posterior from 0.23 / 0.33 / 0.44 to **0.31 / 0.30 / 0.39** (c_mix = 0). The constraint is what pushes coastal below terrestrial. Without it the two pure strategies are tied. The repo does not record whether the 70–73 ka range was chosen after seeing the Nefud date; the author should state this.

## 4. Drowned-site hotspots normalised by shelf area (fixed)

Shelf area above −120 m is computed from the model grid's 20-arcmin cells. Ranked by score per 1,000 km² of shelf:

| Region | Share of hotspot score (old ranking) | Shelf km² above −120 m | Score per 1,000 km² |
|---|---|---|---|
| Makran & Indus delta | 0.132 | 64,208 | 0.207 |
| Bay of Bengal & Myanmar | 0.399 | 335,941 | 0.119 |
| W India shelf (Gujarat–Kerala) | 0.164 | 245,022 | 0.067 |
| Sri Lanka & Palk Strait | 0.025 | 47,587 | 0.053 |
| Gulf of Aden & S Arabia | 0.007 | 30,948 | 0.024 |
| Oman & Persian Gulf | 0.052 | 290,639 | 0.018 |
| Sunda shelf (Malaya, Sumatra, Java) | 0.205 | 1,602,868 | 0.013 |
| Red Sea & Bab-el-Mandeb | 0.006 | 183,811 | 0.004 |
| Sahul shelf | 0.002 | 939,527 | 0.000 |

The **Makran & Indus delta** ranks first per unit of shelf. The Bay of Bengal's first place in the raw ranking is partly area. The Sunda shelf drops far down once its very large area is counted.

## Additional problem found in review

- **Front-speed prior.** The model's `v_max` prior (0.86–5.08 km/yr, 10–90%) is far wider than published estimates: about 2 km/yr (Field & Lahr 2005) and at least 1 km/yr (Mellars 2006); see `inputs/other_sources.json` in the Atlas repo. The ABC posterior shifts up (median 3.1 km/yr), so the short chronology leans on the fast tail of an unsourced prior.

## Still open (needs outside data)

- Replace the hand-drawn biome boxes and humidity pulses with published paleoclimate output (Beyer et al. 2020, Sci Data 7:236, figshare 12293345; Krapp et al. 2021, Sci Data 8:228, osf.io/8n43x; Timmermann & Friedrich 2016 at climatedata.ibs.re.kr).
- Replace the 17-point sea-level curve with the cited reconstruction (Spratt & Lisiecki 2016, NOAA NCEI; Grant et al. 2012, NCEI study 13685).
- Source the front-speed prior (see above).
- Compare against the least-cost-path literature (Field & Lahr 2005; Field et al. 2007; other studies listed in the Atlas `inputs/other_sources.json`, item G).

## Same tests on the LOVECLIM climate (14 November 2026)

Habitat and re-tooling classes from the LOVECLIM transient run behind Timmermann & Friedrich 2016 (`code/model_ibs.py`, runner `code/run_loveclim.sh`, output `results/loveclim/review_ibs.json`). Long chronology, 400 paired draws.

| c_mix | P(mixed first at Sahul) | Posterior coastal / megafauna / mixed | Without the Nefud test |
|---|---|---|---|
| 0 | 0.99 | 0.16 / 0.36 / 0.48 | 0.33 / 0.26 / 0.41 |
| 0.1 | 0.82 | 0.17 / 0.38 / 0.45 | 0.34 / 0.27 / 0.39 |
| 0.2 | 0.19 | 0.18 / 0.41 / 0.41 | 0.35 / 0.29 / 0.36 |

Coastal reaches Sahul before megafauna in 82.5% of runs (75% on the hand-drawn climate). Both review conclusions hold under the modelled climate: mixed leads only when flexibility is free, and the Nefud test decides the coastal vs megafauna ranking. Without it, coastal fits better than megafauna.
