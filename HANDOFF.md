# Handoff: next phase of "Coast or herds?"

**Owner:** Trey (GitHub: TreyModels). **Status:** working model and four experiments done (`FINDINGS.md`). The next phase needs **network access** the original workspace didn't have. Its allowlist blocked climatedata.ibs.re.kr, Zenodo, NOAA and doi.org.

## Goal

Turn the *Homo erectus* walk / raft / canoe result into something publishable. Re-run it on peer-reviewed paleoclimate habitat layers instead of hand-drawn biome boxes, with better coastlines, crossings and site data. Then do the same upgrade for the *sapiens* runs if time allows.

**The claim to test:** sites on Sulawesi (≥1.04 Ma), Flores (1.02 Ma) and Luzon (0.71 Ma), together with the **absence** of any pre-*sapiens* humans in Sahul, the Andamans and Socotra, imply rare accidental crossings: about 1 per 50–170 kyr. They rule out walking only and deliberate canoe voyaging.

## Tasks, in order

### 1. Get the climate–habitat data (highest priority)

- **Dataset:** ICCP/IBS transient CESM1.2 simulation and hominin habitat-suitability layers.
  - Landing pages: https://climatedata.ibs.re.kr/data/hominin-habitats and https://climatedata.ibs.re.kr/data/hominin-habitats/2ma-transient-climate-simulation-data
  - 3 Ma extension: https://climatedata.ibs.re.kr/data/3ma-transient-climate-simulation (needed, because Shangchen is 2.12 Ma).
  - Data DOI 10.22741/iccp.20230001.
  - Papers: Timmermann et al. 2022, *Nature* 604:495 (doi:10.1038/s41586-022-04600-9); Yun et al. 2023, *Clim. Past* 19:1951.
- **Format:** 3.75° grid. Habitat time axis 0001–2001 = 2 Ma → present in 1-kyr steps; NaN where a species is absent. Served via OPeNDAP and LAS. The OPeNDAP URL sits behind the "Data Access" button, so copy it from the page. The page states no licence: cite both papers and email the authors before publishing.
- **Fetch the Eurasian *H. erectus* habitat variable** for lon 20–165, lat −25–50:
  ```bash
  pip install xarray netCDF4 scipy
  python -c "from climate_inputs import fetch_habitat; fetch_habitat('<OPeNDAP URL>', 'erectus_habitat.nc')"
  ```
- **Also pull, if available:** precipitation, temperature and NPP from the 3 Ma run, for the >2 Ma part and for *sapiens*.
- **Run with it:**
  ```bash
  cd code
  ERECTUS_HABITAT=../data/erectus_habitat.nc ERECTUS_HABITAT_VAR=<var name> python erectus.py 400
  python er_analyze.py && python er_map.py
  ```
  The hook already exists and has been tested on synthetic data (`climate_inputs.py selftest`). With the variable set, terrestrial habitat = simulated suitability × altitude penalty; sea level, crossings and scoring are unchanged.
- **Caution:** at 3.75° (~400 km) the habitat field is far coarser than the 20′ landscape. Keep coastlines from bathymetry, and treat the habitat as a smooth climate field. Check that islands don't inherit ocean NaNs (zero suitability); fill coastal cells from the nearest land value if they do.

### 2. Better coastlines

- Replace ETOPO 10′ (via pyvista/vtk-data) with **GEBCO 2024** or **ETOPO 2022 at 1′–30″**, resampled to ~5′.
- Use a GIA-aware relative sea level where possible (Dobson/Spada/Galassi 2024; Lambeck regional curves).
- Use the **Clark et al. 2025** global sea-level curve (Science 390, eadv8389) as the main history. The classic curve becomes a sensitivity check.
- Re-check the key straits at each lowstand: Bab-el-Mandeb, Hormuz, Lombok, Makassar, Timor–Sahul, Korea/Tsushima, Mindoro/Palawan–Luzon.

### 3. Crossings

- Replace the `w0·exp(width/λ)` rule with drift-based crossing probabilities, following **Bird et al. 2019**, *Sci. Rep.* 9:8220 (doi:10.1038/s41598-019-42946-9; 17 Wallacean crossing points, ocean drift and visibility).
- Ask the Bradshaw/Bird group whether their drift matrices can be shared.
- Keep three modes: walk, accidental raft (passive drift only), deliberate craft (minimal headway, e.g. 0.5 kn).

### 4. Site database

- Replace the hand-picked 10 sites with a systematic pull from **ROAD** (ROCEEH Out of Africa Database, https://www.roceeh.uni-tuebingen.de/roadweb) or the Timmermann et al. species database (supplement).
- Keep explicit age uncertainty: sample site ages, don't use point minima.
- **Absence list:** Sahul, the Andamans, Socotra, Sri Lanka (check), Japan (flag as contested), Ryukyus, Taiwan (Penghu jaw = Denisovan; check timing).
- Score with a proper likelihood. Each site gets P(reached by its sampled age). Each absence gets a detection-adjusted probability, i.e. how likely we'd have found sites if people were there.

### 5. Inference and robustness

- Approximate Bayesian computation over crossing rate, speed and biome cost; report posteriors, not just hit rates.
- Leave-one-site-out validation and cross-validation on the island sites.
- Sensitivity checks: grid resolution (20′ vs 5′), sea-level history, habitat source (boxes vs CESM).
- Write down the predictions (e.g. Timor or Sumba empty or occupied before 0.5 Ma) before looking at any new data.

### 6. *Sapiens* upgrade (optional)

- Same climate swap, using the 3 Ma run's last 100 kyr, or the Krapp et al. 2021 / Beyer et al. 2021 high-resolution emulators for 120 ka → present.
- Have a population geneticist redo the genetics with a real spatial coalescent simulator (e.g. msprime with a stepping-stone deme graph built from the colonisation tree). The `genetics.py` path-sum is only an approximation.

## Deliverables

1. `data/` with the fetched layers and a `DATA_SOURCES.md` (URLs, access date, citation, licence status).
2. Re-run results: `results/er_cesm.json`, an updated `FINDINGS.md`, updated figures and page.
3. A 3,000–5,000-word preprint draft for SocArXiv or OSF. Target journals: *J. Human Evolution*, *Quaternary Science Reviews*, *Sci. Reports*.
4. A short pitch email to a domain group: Bradshaw/Bird (Sahul and crossings), Timmermann (climate–habitat) or Brumm (Wallacea sites).

## Done means

- The *erectus* result is re-run on CESM habitat with ≥400 draws and reported both ways (boxes vs CESM).
- Each verdict holds or fails: walk = 0, canoe overshoots, rare rafting fits. Report it either way.
- Every number in the preprint traces to a script plus a results file in the repo.

## Where things are

- **Engine:** `code/model.py` (landscape, sea level, crossings, time-dependent Dijkstra in numba).
- **Sapiens:**
  - `code/mc.py`: Monte Carlo.
  - `code/figs.py`, `code/frames.py`, `code/climate_frames.py`: figures and slider.
  - `code/exp.py`: evidence fitting (ABC) and drowned-site experiments.
  - `code/genetics.py`, `code/gen_analyze.py`: genetics layer.
- **Erectus:** `code/erectus.py`, `code/er_analyze.py`, `code/er_map.py`.
- **New:** `code/climate_inputs.py` (fetch, regrid and hook for CESM habitat).
- **Page:** `coast-or-herds.html` (self-contained); template `code/page_template.html`.
- **Paths:** scripts assume `/home/claude/work/`; change `cross.npz`, `grid.npz` and similar paths to your folder. Regenerate `cross.npz` with `build.py` (~6 s).
