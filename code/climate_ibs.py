"""
LOVECLIM climate forcing of Timmermann & Friedrich (2016) on the model's 20' grid.

Source: data/external/ibs_hds/climate-forcing.nc (IBS ICCP; 1° global, mtime = 1..130 kyr BP)
  netpp     net primary production (kgC m-2 yr-1)
  desert_fr desert fraction (%)
  temp      surface air temperature (°C)

Ocean / missing cells are filled from the nearest valid land cell before regridding, so shelf that is
dry at glacial sea level inherits the climate of the adjacent coast (LOVECLIM's land mask is coarse).
Slices KA_LO..KA_HI are kept so the Monte Carlo can shift the climate chronology by a few kyr.

Run once to build the cache:  python climate_ibs.py
"""
import os, numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
import model as M

SRC = os.path.join(M.ROOT, 'data', 'external', 'ibs_hds', 'climate-forcing.nc')
CACHE = os.path.join(M.ROOT, 'data', 'external', 'ibs_hds', 'climate_20min.npz')
KA_LO, KA_HI = 30, 105                      # kyr BP kept in the cache (model window is 100-35)
KA = np.arange(KA_LO, KA_HI + 1)            # ascending
VARS = ('netpp', 'desert_fr', 'temp')
# Unresolved-land screen: LOVECLIM cells that are almost desert-free yet almost barren are physically
# inconsistent. They are small islands / thin coasts that the coarse model barely resolves (ocean-diluted),
# e.g. Nusa Tenggara. They are treated as missing and filled from the nearest resolved land, like ocean cells.
# Only WARM cells are screened: cold steppe / tundra / high plateau really is low-NPP without being desert.
UNRES_NPP, UNRES_DES, UNRES_T = 0.15, 30.0, 18.0   # netpp < 0.15 AND desert_fr < 30 % AND temp >= 18 °C


def _fill_nearest(a):
    bad = ~np.isfinite(a)
    if not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def build():
    from netCDF4 import Dataset
    ds = Dataset(SRC)
    lat = np.asarray(ds['lat'][:], float); lon = np.asarray(ds['lon'][:], float)
    mt = np.asarray(ds['mtime'][:], float)
    # region with a margin so interpolation and nearest-fill near the edges use real data
    la_sel = (lat >= M.LAT.min() - 6) & (lat <= M.LAT.max() + 6)
    lo_sel = (lon >= M.LON.min() - 6) & (lon <= M.LON.max() + 6)
    la, lo = lat[la_sel], lon[lo_sel]
    o_la = np.argsort(la); la = la[o_la]
    pts = np.column_stack([M.LAT2.ravel(), M.LON2.ravel()])

    def raw(v, ti):
        a = np.ma.filled(np.ma.masked_invalid(ds[v][ti][np.ix_(la_sel, lo_sel)].astype(float)), np.nan)
        a = np.where(np.abs(a) > 1e20, np.nan, a)
        if v == 'temp':
            a = np.where(a <= -99, np.nan, a)
        return a[o_la]

    out = {v: np.zeros((len(KA), M.NY, M.NX), np.float32) for v in VARS}
    unres = np.zeros((len(KA), len(la), len(lo)), bool)
    for k, ka in enumerate(KA):
        ti = int(np.argmin(np.abs(mt - ka)))
        assert abs(mt[ti] - ka) < 1e-6, f'no slice for {ka} ka'
        R = {v: raw(v, ti) for v in VARS}
        bad = (R['netpp'] < UNRES_NPP) & (R['desert_fr'] < UNRES_DES) & (R['temp'] >= UNRES_T)
        unres[k] = bad
        for v in VARS:
            a = np.where(bad, np.nan, R[v])
            a = _fill_nearest(a)
            f = RegularGridInterpolator((la, lo), a, method='linear', bounds_error=False, fill_value=None)
            out[v][k] = f(pts).reshape(M.NY, M.NX)
    out['unres_src'] = unres
    out['src_lat'] = la; out['src_lon'] = lo
    land = ~M.water_mask(-60.0)                                     # typical glacial land in the window
    win = (KA >= 35) & (KA <= 100)
    npp_land = out['netpp'][win][:, land]
    out['NPP90'] = np.float32(np.percentile(npp_land, 90))         # reference 'rich savanna' productivity
    out['ka'] = KA
    np.savez_compressed(CACHE, **out)
    return out


_C = None
def load():
    global _C
    if _C is None:
        if not os.path.exists(CACHE):
            build()
        z = np.load(CACHE)
        _C = {k: z[k] for k in z.files}
    return _C


def at(bp, shift_yr=0.0):
    """(npp, desert_fr %, temp °C) 2-D fields at years BP, chronology shifted by shift_yr (linear in time)."""
    C = load()
    ka = np.clip((bp + shift_yr) / 1000.0, KA_LO, KA_HI)
    x = ka - KA_LO; i0 = int(np.floor(x)); i1 = min(i0 + 1, len(KA) - 1); w = x - i0
    return tuple((1 - w) * C[v][i0] + w * C[v][i1] for v in VARS)


if __name__ == '__main__':
    C = build()
    land = ~M.water_mask(-60.0)
    print('cache:', CACHE)
    print('NPP90 (kgC/m2/yr, land, 100-35 ka):', float(C['NPP90']))
    for v in VARS:
        a = C[v][(C['ka'] >= 35) & (C['ka'] <= 100)][:, land]
        print(f'{v:10s} land p5/p50/p95: {np.percentile(a, 5):.3f} {np.percentile(a, 50):.3f} {np.percentile(a, 95):.3f}')
    # Saharo-Arabian belt: does desert fraction pulse with the known humid phases?
    belt = land & (M.LAT2 > 16) & (M.LAT2 < 30) & (M.LON2 > 35) & (M.LON2 < 58)
    print('Arabia mean desert_fr % by ka:')
    print(' '.join(f'{int(k)}:{C["desert_fr"][i][belt].mean():.0f}' for i, k in enumerate(C['ka']) if k % 5 == 0))
    # where the unresolved-land screen fired (source 1° cells, any time 35-100 ka)
    U = C['unres_src'][(C['ka'] >= 35) & (C['ka'] <= 100)].any(0)
    sla, slo = np.meshgrid(C['src_lat'], C['src_lon'], indexing='ij')
    print('unresolved-land source cells (any time 35-100 ka):', int(U.sum()))
    regions = {'Nusa Tenggara/Timor': (114, 127, -11, -7), 'Sulawesi/Moluccas': (118, 131, -7, 3),
               'Sunda (Sumatra/Java/Borneo)': (95, 119, -9, 7), 'New Guinea/N Sahul': (127, 150, -15, 0),
               'Arabia/Red Sea': (32, 60, 10, 31), 'India/Sri Lanka': (66, 92, 5, 26)}
    for name, (a, b, c, d) in regions.items():
        m = U & (slo >= a) & (slo <= b) & (sla >= c) & (sla <= d)
        print(f'  {name:28s} {int(m.sum())}')
    # Nusa Tenggara check: big-game NPP after the fix
    nt = (M.LAT2 > -10.5) & (M.LAT2 < -7.5) & (M.LON2 > 115) & (M.LON2 < 126) & ~M.water_mask(-70.0)
    print('Nusa Tenggara land netpp after fix, p10/p50/p90 over 35-100 ka:',
          np.round(np.percentile(C['netpp'][(C['ka'] >= 35) & (C['ka'] <= 100)][:, nt], [10, 50, 90]), 3))
