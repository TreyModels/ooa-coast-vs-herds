# Findings so far

All numbers are Monte Carlo shares of runs, with paired parameter draws across strategies (same climate, sea level, speed and noise for each strategy). Engine: time-dependent least-cost dispersal on a 20′ (~37 km) ETOPO paleo-landscape. See `README.md` for the code map, and `coast-or-herds.html` for charts and maps.

Strategies:
- **Coastal:** shore-bound, with fresh water from coastal springs in dry phases.
- **Megafauna / terrestrial:** follows grassland and game, pays a re-tooling delay on entering each new biome.
- **Mixed:** takes whichever is better in each cell.

## 1. *Homo sapiens* out of Africa (exit 85–72 ka, 400 paired draws)

| | Coastal | Megafauna | Mixed |
|---|---|---|---|
| Horn → Sahul transit, median (kyr) | 12.2 | 13.3 | 10.9 |
| Share of route now under the sea (median) | 46% | 33% | 39% |
| Reaches Sulawesi by 67.8 ka (Muna rock art) | 50% | 48% | 62% |
| Reaches Sahul by 65 ka (Madjedbebe) | 56% | 51% | 62% |
| Occupies interior Arabia 85–70 ka (Nefud) | 59% | 95% | 96% |
| Sulawesi in time **and** interior Arabia | 33% | 48% | 62% |
| Short chronology: Near-East start 51–47 ka, Sahul by 45 ka | 0.7% | 0.5% | 5% |

- **Head-to-head:** coastal reaches Sahul before megafauna in **75%** of paired runs, by a median of ~1,100 years.
- **Decisive assumption:** coastal productivity relative to good savanna. Coastal wins 36% of runs in the lowest quartile (0.70–0.84×) and 95% in the highest (1.08–1.20×). It wins 98% when the coast is rich, keeps dry-phase water and biome changes are costly, and 23% in the opposite corner.
- **Evidence fit (rejection ABC on Sulawesi + Sumatra + Nefud):**
  - Posterior by strategy: coastal 23%, megafauna 33%, mixed 44%.
  - Fitting needs front speed ≳ 2 km/yr (median ~3 km/yr) and an early exit (79.5–84.5 ka fits 80–84% of runs vs 25% for 72–74.5 ka).
  - The evidence does **not** constrain coastal productivity: posterior = prior (median 0.98×). The current dry-land record can't settle coast vs herds.
- **Where the drowned route lies (coastal + mixed):**
  - Share of the drowned route by region: Bay of Bengal–Myanmar shelf 40%, Sunda shelf 21%, W India shelf 16%, Makran–Indus 13%.
  - Median present depth ~24 m (10–90%: 5–52 m) at cell scale.
  - The Red Sea and S Arabian shelves are below grid resolution.

## 2. Genetics layer (main wave from Persian-Plateau hub, 50–46 ka start, 150 draws)

| | Coastal | Megafauna | Mixed |
|---|---|---|---|
| Papuans sister to East Asians, before S India | **100%** | 55% | 94% |
| Spread of East Eurasian split times, median (kyr) | 2.4 | 3.7 | 1.9 |

- **Branching order:** every strategy gives the east/west split first (~48 ka) with near-star eastern splits (~46–47 ka). Only coast-anchored spread reliably reproduces the observed Papuan–East Asian clade (Lipson & Reich 2017).
- **Serial-founder diversity loss:** all strategies give the observed west→east decline, lowest in Sahul. The coastal 1-D front roughly doubles drift, but the absolute scale depends on deme size, so this test doesn't discriminate.
- **Ghost of an earlier wave:** if an early wave (e.g. the ≥67.8 ka Sulawesi painters) still occupied the route, Papuan ancestors cross 211–363 cells (~7,800–13,400 km) of its land. Keeping early-wave ancestry ≤2% (Pagani et al. 2016's upper estimate; others find ~0) needs < ~1 interbreeding per 10,000–18,000 per 37-km cell. ⇒ The early wave was gone or nearly gone by ~50 ka, or wasn't *H. sapiens*. This holds for any strategy combination.

## 3. *Homo erectus*-grade dispersal (2.6 → 0.1 Ma, 400 draws, 3 strategies × 3 crossing modes, 2 sea-level histories)

Fit = share of runs reaching all 10 dated sites on time **and** staying out of Sahul, the Andamans and Socotra. Japan is excluded: a glacial land bridge is debated.

| Crossing ability | Coastal | Megafauna | Mixed |
|---|---|---|---|
| Walk only | 0% | 0% | 0% |
| Rare accidental rafting | 2% | **10%** | **10%** |
| Dugout canoes | 0% | 0% | 0% |

- **Not beachcombers:**
  - Shangchen by 2.12 Ma: coastal 71%, megafauna 99%, mixed 100%.
  - Dmanisi by 1.85 Ma: coastal 37% vs ~100% for the other two.
  - All seven mainland sites on time: coastal 33% vs 99–100%.
- **Walking can't explain the islands:** Sulawesi (≥1.04 Ma), Flores (1.02 Ma) and Luzon (0.71 Ma) are reached in 0% of walk-only runs.
- **Canoes overshoot:** they reach the islands (86%), but also reach Sahul in 86% of runs (usually almost immediately), the Andamans in 82% and Socotra in ~46%. None of these has pre-*sapiens* evidence.
- **Rare rafting fits:** the best rate is ~1 successful crossing per 50–170 kyr; 29% of runs at 100–200 kyr fit everything.
- **Sea-level history (classic vs Clark et al. 2025)** barely changes any result.

## Main caveats

- Biomes are hand-drawn boxes, not paleoclimate output. This is the #1 thing to replace; see `HANDOFF.md`.
- Sea-crossing widths are approximated on the 20′ grid (≈25 km correction); no drift or current modelling.
- There's no explicit demography (population size, Allee effects, competition with archaics). Genetics uses path-based approximations, not coalescent simulation.
- Present depth is a cell median, so narrow drowned beaches are under-counted.
- Results show the logic of each strategy; they don't date any event.
