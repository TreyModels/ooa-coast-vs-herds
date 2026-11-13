# Coast or herds? Out-of-Africa dispersal model

**Start here:** `FINDINGS.md` (what we found) · `HANDOFF.md` (next-phase assignment, incl. fetching the Timmermann/IBS climate–habitat data) · `coast-or-herds.html` (full write-up with maps and charts).

Time-dependent least-cost dispersal model comparing three foraging strategies
(coastal, megafauna/terrestrial, mixed) from the Horn of Africa to Sahul, 100–35 ka.

## Contents
- `coast-or-herds.html` – the full write-up with maps, climate time slider and charts (open in any browser, works offline apart from fonts).
- `code/` – the model and analysis scripts.
- `data/grid.npz` – 20′ study grid (lon 25–160E, lat 20S–45N) derived from ETOPO 10′ ice-surface relief (zmin / zmed / zmax per cell).
- `results/monte_carlo_runs.csv` – one row per run (400 long-chronology "L" + 400 short-chronology "S"): sampled parameters plus arrival time (years BP; -inf = never) at each site and at Sahul for each strategy, route-sunk fraction, number of sea crossings.
- `results/density.npz` – how many runs used each cell on the fastest route to Sahul, and occupancy at the snapshot date.
- `results/summary.json`, `results/exp.json` – the numbers shown on the page.
- `figures/` – static maps.

## Re-running
Python 3 with numpy, scipy, pandas, matplotlib, numba, pillow (and pyvista only to rebuild the grid).
Run every script from an empty working folder next to `code/`, for example:

    mkdir work && cd work
    python ../code/build.py
    python ../code/mc.py 400

Outputs (`cross.npz`, `mc.pkl`, `density.npz`, figures) are written to that folder; the grid is read from `data/`.

1. `prep.py` – builds `grid.npz` from ETOPO (needs `etopo10.npy`; regenerate from pyvista/vtk-data `Data/EarthModels/ETOPO_10min_Ice.vtp`). Skip it and use the included grid.
2. `build.py` – precomputes sea-crossing edges per sea-level bin → `cross.npz` (~6 s).
3. `mc.py 400` – paired Monte Carlo, both chronologies (~6 min on 2 cores) → `mc.pkl`, `density.npz`.
4. `frames.py` – 150 runs keeping full arrival grids (for the climate slider) → `arr_L.npy` (~140 MB).
5. `figs.py`, `climate_frames.py`, `exp.py` – figures, slider frames, ABC and drowned-site experiments.
6. `analyze.py`, `analyze2.py` – text summaries and sensitivity tables.

## Key assumptions (all sampled with uncertainty in `model.sample_params`)
Sea-level curve (approximate composite, ±8 m, ±1.5 kyr); Saharo-Arabian humidity pulses (MIS 5c, 5a, ~55 ka);
hand-drawn biome zones with river corridors and the exposed Gulf basin; front speed 0.7–6 km/yr;
coastal productivity 0.7–1.2× best savanna; dry-phase coastal fresh water; biome re-tooling cost 0–1,500 yr;
sea crossings with waiting time w0·exp(width/λ). This is a model of the logic of each strategy, not a date for any event.

## Genetics layer
- `genetics.py` – 150 runs per strategy of the main wave (ancestors of living non-Africans) spreading from the Persian-Plateau hub, 50–46 ka start, to 25 ka (it sets `OOA_T_END=25000` itself). Saves colonisation paths and front widths to `wave2.pkl`.
- `gen_analyze.py` – split times and tree topology between reference populations, serial-founder drift (Σ 1/front width along the founding path), and ghost ancestry picked up from an earlier wave (`arr_L.npy`, from `frames.py`) → `results/gen.json`.

## Homo erectus layer
- `erectus.py N` – 2.6 → 0.1 Ma on the same engine (5-kyr steps), 3 strategies × 3 sea-crossing abilities (walk / rare rafting / dugout canoes), two sea-level histories (classic vs Clark et al. 2025). Scored against 10 dated sites and 4 never-reached areas → `erectus.pkl`, `erectus_reach.npz`.
- `er_analyze.py` → `results/er.json`; `er_map.py` → `figures/map_erectus.png`; per-run results in `results/erectus_runs.csv`.
- `climate_inputs.py` – fetch + regrid hook for the Timmermann et al. 2022 CESM habitat layers. Set `ERECTUS_HABITAT=path.nc` (and `ERECTUS_HABITAT_VAR`) to run `erectus.py` on them instead of the hand-drawn biomes. `python ../code/climate_inputs.py` runs a self-test on synthetic data.

Full run order from an empty `work/` folder: `build.py` → `mc.py 400` → `figs.py` → `frames.py` → `climate_frames.py` → `exp.py` → `genetics.py` → `gen_analyze.py` → `erectus.py 400` → `er_analyze.py` → `er_map.py`.

## Attribution and reuse

Credit and cite **Dr. Trey Menefee** when using, reproducing, sharing, or adapting his material. Noncommercial use is permitted under the conditions in [LICENSE](LICENSE); commercial use requires prior permission. Third-party materials retain their own terms.
