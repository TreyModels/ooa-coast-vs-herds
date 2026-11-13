"""
Plug real paleoclimate output into the erectus / sapiens runs in place of the hand-drawn biome boxes.

Source: ICCP/IBS 2-3 Myr transient CESM1.2 simulation and its hominin habitat-suitability layers
  Timmermann et al. 2022, Nature 604:495-501, doi:10.1038/s41586-022-04600-9
  Yun et al. 2023, Clim. Past 19:1951, doi:10.5194/cp-19-1951-2023 (data doi 10.22741/iccp.20230001)
  Landing pages:
    https://climatedata.ibs.re.kr/data/hominin-habitats
    https://climatedata.ibs.re.kr/data/hominin-habitats/2ma-transient-climate-simulation-data
    https://climatedata.ibs.re.kr/data/3ma-transient-climate-simulation
  Grid 3.75 deg x 3.75 deg; habitat time axis 0001..2001 = 2 Ma .. present in 1-kyr steps; NaN where a
  species did not exist.

Usage
  1) fetch_habitat(url, out)      # OPeNDAP URL copied from the "Data Access" button on the landing page
  2) H = HabitatField(out, var)   # var = the Eurasian H. erectus suitability variable in that file
  3) erectus.py: set ERECTUS_HABITAT=/path/to/file.nc (and ERECTUS_HABITAT_VAR) and the run uses
     H.at(t_bp) as terrestrial habitat quality instead of the biome boxes. Everything else is unchanged.
"""
import numpy as np


def fetch_habitat(opendap_url, out_path, lon=(20, 165), lat=(-25, 50)):
    """Subset the study window over OPeNDAP and save locally as netCDF (needs network access)."""
    import xarray as xr
    ds = xr.open_dataset(opendap_url)
    lonname = 'lon' if 'lon' in ds.coords else [c for c in ds.coords if c.lower().startswith('lon')][0]
    latname = 'lat' if 'lat' in ds.coords else [c for c in ds.coords if c.lower().startswith('lat')][0]
    lat_sorted_up = float(ds[latname][0]) < float(ds[latname][-1])
    sub = ds.sel({lonname: slice(*lon), latname: slice(*(lat if lat_sorted_up else lat[::-1]))})
    sub.load().to_netcdf(out_path)
    print('saved', out_path, dict(sub.sizes), list(sub.data_vars))
    return out_path


class HabitatField:
    """Suitability (0..1) on the model's 20' grid at any time BP; linear in time, bilinear in space."""

    def __init__(self, path, var, time_origin_bp=2_000_000, step_years=1000):
        import xarray as xr
        import model as M
        ds = xr.open_dataset(path, decode_times=False)
        da = ds[var]
        tname = [d for d in da.dims if d.lower().startswith('t')][0]
        latn = [d for d in da.dims if d.lower().startswith('lat')][0]
        lonn = [d for d in da.dims if d.lower().startswith('lon')][0]
        da = da.transpose(tname, latn, lonn)
        self.t_bp = time_origin_bp - (np.arange(da.sizes[tname])) * step_years   # 0001 -> 2 Ma
        self.lat = da[latn].values; self.lon = da[lonn].values
        self.data = np.nan_to_num(da.values.astype(np.float32), nan=0.0)
        v = self.data[self.data > 0]
        if v.size and v.max() > 1.0:                                  # rescale to 0..1 if needed
            self.data = self.data / np.nanpercentile(v, 99)
        self.M = M
        from scipy.interpolate import RegularGridInterpolator
        self._RGI = RegularGridInterpolator

    def at(self, t_bp):
        M = self.M
        i = np.clip(np.searchsorted(-self.t_bp, -t_bp), 1, len(self.t_bp) - 1)
        t0, t1 = self.t_bp[i - 1], self.t_bp[i]
        w = 0.0 if t0 == t1 else (t0 - t_bp) / (t0 - t1)
        f = (1 - w) * self.data[i - 1] + w * self.data[i]
        lat, lon, ff = self.lat, self.lon, f
        if lat[0] > lat[-1]: lat, ff = lat[::-1], ff[::-1]
        rgi = self._RGI((lat, lon), ff, bounds_error=False, fill_value=0.0)
        pts = np.c_[M.LAT2.ravel(), M.LON2.ravel()]
        return np.clip(rgi(pts).reshape(M.NY, M.NX), 0, 1)


def selftest():
    """Synthetic file -> HabitatField -> 20' grid, to check the plumbing without network."""
    import xarray as xr, tempfile, os
    lat = np.arange(-24.375, 50, 3.75); lon = np.arange(20.625, 165, 3.75); t = np.arange(1, 2002)
    data = np.random.default_rng(0).random((len(t), len(lat), len(lon))).astype(np.float32)
    ds = xr.Dataset({'erectus': (('time', 'lat', 'lon'), data)}, coords={'time': t, 'lat': lat, 'lon': lon})
    p = os.path.join(tempfile.gettempdir(), 'habitat_selftest.nc'); ds.to_netcdf(p)
    H = HabitatField(p, 'erectus'); a = H.at(1_500_000)
    print('selftest ok', a.shape, float(a.min()), float(a.max()))


if __name__ == '__main__':
    selftest()
