# IBS Human Dispersal Simulation (125 ka) — local backup

Backup copy of the model output from:

> Timmermann, A. & Friedrich, T. (2016). Late Pleistocene climate drivers of early human migration.
> *Nature* 538, 92–95. https://doi.org/10.1038/nature19365

Downloaded 2026-11-13 from the IBS Center for Climate Physics (ICCP) climate data server:

- Dataset page: https://climatedata.ibs.re.kr/data/human-dispersal-simulation
- Files: https://climatedata.ibs.re.kr:8443/thredds/catalog/human-dispersal-simulation/catalog.html
  (THREDDS `fileServer/human-dispersal-simulation/<file>`; also OPeNDAP at `:9090/dods/public-data/human-dispersal-simulation`)

As of the download date this is the only public host. There is no DOI and no copy on Zenodo, PANGAEA,
NOAA NCEI, Dryad or figshare. The one partial re-host is `population_data/timmermann.npz` (density only,
NumPy format) in github.com/rwalker1501/Cultural-Epidemiology.

## Files

All on a global 1° grid (lat 180 × lon 360; lon −19.5…339.5).

| File | Variables | Notes |
|---|---|---|
| `arrival-time.nc` | `arrivaltime_a`…`_f` (ka BP) | Time since last continuous settlement. a = early exit, b = late exit, c = early exit without D-O variability, d/e/f = constant climate at 105/70/21 ka |
| `climate-forcing.nc` | `temp` (°C), `netpp` (kgC m⁻² yr⁻¹), `desert_fr` (%) | 130 × 1-kyr steps; LOVECLIM transient run plus millennial variability; the forcing behind the dispersal model |
| `human_density.nc` | `human_density_a`, `human_density_b` (individuals per 100 km²) | 1-kyr snapshots, scenario A (early exit) and B |

The server also lists `human-density.nc`, `human-density_old.nc` and `human-density-old.nc`, all the same size as
`human_density.nc`. They were not downloaded.

Checksums: `SHA256SUMS` (verify with `shasum -a 256 -c SHA256SUMS`).

## License

The ICCP site states that data on its server are shared under "Creative Commons Attribution 4.0 International
License, CC BY-NC", which contradicts itself: CC BY 4.0 allows commercial use and CC BY-NC does not. Treat it as
**attribution required, noncommercial use only**. Cite the dataset page and doi:10.1038/nature19365. Ask ICCP
(dataset contact: Elke Zeller) before redistributing these files or using them commercially.

These files are not part of this project's own license and should not be committed to the GitHub repo:
`human_density.nc` (130 MB) is over GitHub's 100 MB per-file limit.
