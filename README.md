# Coast or herds? Out-of-Africa dispersal model

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
