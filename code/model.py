"""
Out-of-Africa dispersal: coastal 'beachcomber' vs terrestrial 'megafauna-tracking' foragers.

Time-dependent least-cost (earliest-arrival) model on a 20' paleo-landscape:
  * land/sea from ETOPO relief + a time-varying sea-level curve
  * Saharo-Arabian humidity pulses ('Green Arabia') and monsoon strength vary through time
  * each strategy turns the landscape into a front-speed surface (km/yr); cells with speed 0
    are uninhabitable at that time and a population cannot wait in them
  * sea gaps are crossed with an expected waiting time that grows with gap width
Monte Carlo draws every uncertain parameter; strategies are compared on identical draws.
"""
import numpy as np, heapq
from scipy import ndimage
import numba as nb

import os
G = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'grid.npz'))
ZMAX, ZMED, ZMIN, LAT, LON = G['zmax'], G['zmed'], G['zmin'], G['lat'], G['lon']
NY, NX = ZMAX.shape
N = NY * NX
LAT2, LON2 = np.meshgrid(LAT, LON, indexing='ij')
CELL_KM = 111.32 / 3

# ---------------------------------------------------------------- time axis
T_START, T_END, DT = 100_000, 35_000, 1000          # years BP
K = (T_START - T_END) // DT + 1
T_K = T_START - np.arange(K) * DT                    # BP at start of slice k

def bp_to_tau(bp):  return T_START - bp
def tau_to_bp(tau): return T_START - tau

# ---------------------------------------------------------------- sea level (m), approximate composite
# Shape follows Red Sea / global RSL reconstructions (Grant et al. 2014; Spratt & Lisiecki 2016);
# MC adds an offset and a time jitter.
SL_KA = np.array([100, 95, 90, 85, 80, 76, 72, 68, 64, 60, 57, 54, 51, 48, 45, 40, 35])
SL_M  = np.array([-25, -40, -50, -30, -25, -40, -60, -80, -88, -80, -65, -72, -78, -75, -85, -80, -85])

def sea_level(bp, offset=0.0, jitter=0.0):
    return np.interp((bp + jitter) / 1000, SL_KA[::-1], SL_M[::-1]) + offset

# ---------------------------------------------------------------- static landscape
def _lakes_and_depressions():
    present_water = (ZMED < 0) | (ZMIN < -30)
    lab, n = ndimage.label(present_water, structure=np.ones((3, 3)))
    border = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])) - {0}
    ocean_now = np.isin(lab, list(border))
    caspian = (LAT2 > 36.4) & (LON2 > 46.5) & (LON2 < 55.5) & (ZMED < -27)
    black = (LAT2 > 40.8) & (LON2 > 27) & (LON2 < 42) & (ZMED < 0)
    gulf = (LAT2 > 23) & (LAT2 < 31) & (LON2 > 47) & (LON2 < 57.5) & (ZMED < 0)   # Persian Gulf is sea, not a basin
    depression = present_water & ~ocean_now & ~caspian & ~gulf       # Qattara, Danakil, Dead Sea area...
    return depression, caspian | black
DEPRESSION, LAKE = _lakes_and_depressions()

def water_mask(sl):
    # 20' cell is sea if most of it is below sea level, or if it holds a channel >30 m deeper than sea level
    w = (((ZMED < sl) | (ZMIN < sl - 30)) & ~DEPRESSION) | LAKE
    # keep small islands (any 10' sub-cell above sea level) when they are not touching a larger landmass,
    # so stepping-stone islands in Wallacea survive the coarse grid without closing real straits
    island = w & ~LAKE & (ZMAX >= sl) & ~ndimage.binary_dilation(~w, structure=np.ones((3, 3)))
    return w & ~island

def _polyline_mask(pts, width=1):
    m = np.zeros((NY, NX), bool)
    for (la0, lo0), (la1, lo1) in zip(pts[:-1], pts[1:]):
        for s in np.linspace(0, 1, 200):
            la, lo = la0 + s * (la1 - la0), lo0 + s * (lo1 - lo0)
            i = int(round((LAT[0] - la) * 3)); j = int(round((lo - LON[0]) * 3))
            if 0 <= i < NY and 0 <= j < NX: m[i, j] = True
    return ndimage.binary_dilation(m, iterations=width)

RIVERS = {
    'Nile':   [(15.6, 32.5), (17.5, 33.9), (19.2, 30.5), (21.8, 31.3), (24.1, 32.9), (25.7, 32.6), (27.2, 31.2), (30.0, 31.2), (31.4, 31.0)],
    'Euphrates': [(37.0, 38.0), (36.0, 38.3), (35.0, 40.4), (34.0, 42.5), (33.0, 44.0), (31.0, 46.0), (30.0, 48.0)],
    'Tigris': [(37.5, 41.0), (36.3, 43.1), (34.0, 44.5), (32.5, 45.5), (31.0, 47.4), (30.0, 48.0)],
    'Gulf paleo-river (Ur-Schatt)': [(30.0, 48.0), (28.8, 49.5), (27.5, 51.5), (26.5, 53.5), (26.2, 55.5), (26.4, 56.5)],
    'Indus':  [(33.5, 72.0), (31.5, 71.0), (29.0, 70.7), (27.0, 68.5), (25.0, 68.2), (24.0, 67.3)],
}
RIVER = np.zeros((NY, NX), bool)
for pts in RIVERS.values(): RIVER |= _polyline_mask(pts, 1)

def _box(lo0, lo1, la0, la1): return (LON2 >= lo0) & (LON2 <= lo1) & (LAT2 >= la0) & (LAT2 <= la1)

# biome classes: 0 hyper-arid belt, 1 semi-arid, 2 moist woodland/savanna, 3 tropical forest, 4 Sahul savanna
BIOME = np.full((NY, NX), 2, np.int8)
BIOME[_box(41, 80, 36, 45) | _box(44, 70, 25, 40) | _box(35.5, 44, 30, 37.5) | _box(25, 38, 10, 16)
      | _box(66, 75.5, 22.5, 32) | (_box(39.5, 52, 0, 12.5) & (ZMED < 1300))] = 1
BIOME[_box(25, 37.5, 16, 31.5) | (_box(34.5, 60, 12, 32.5) & ~_box(25, 36, 31, 33) & ~_box(44, 60, 27.5, 32.5))
      | _box(56, 61, 28, 33) | _box(52, 57, 33, 35.5) | _box(63, 67, 29.5, 31.5) | _box(68.5, 75, 24, 30)] = 0
BIOME[_box(41.8, 46.5, 12.5, 20.5) & (ZMED > 1000)] = 1                # Yemen / Asir highlands
BIOME[_box(93, 158, -10.5, 8.5) | _box(97, 110, 8.5, 14)] = 3           # Sunda, Wallacea, New Guinea
BIOME[(LAT2 < -10.5) & (LON2 > 110)] = 4                                # northern Sahul savanna
BIOME[(LAT2 < -10.5) & (LON2 > 120) & (LON2 < 142) & (LAT2 < -17)] = 1
GULF_BASIN = _box(47, 57, 24, 30.6) & (ZMAX < 0)                        # 'Gulf Oasis' (exposed shelf)
ARID_ZONE = BIOME <= 1

ALT_PEN = np.where(ZMED > 3000, 0.05, np.where(ZMED > 2000, 0.4, np.where(ZMED > 1300, 0.75, 1.0)))

# ---------------------------------------------------------------- neighbours / distances
OFFS = np.array([(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)])
def _nbr_dist():
    d = np.zeros((NY, 8))
    for r in range(NY):
        cl = np.cos(np.radians(LAT[r]))
        for e, (dy, dx) in enumerate(OFFS):
            la = LAT[r] - dy / 3 / 2   # midpoint latitude
            d[r, e] = CELL_KM * np.hypot(dy, dx * np.cos(np.radians(la)))
    return d
NBR_D = _nbr_dist()

# ---------------------------------------------------------------- sea-crossing edges, per sea-level bin
SL_BINS = np.arange(-120, -9, 5)            # -120 ... -10 m
def _crossing_edges(water, max_km=130):
    land = ~water
    coast = land & ndimage.binary_dilation(water, structure=np.ones((3, 3)))
    ci, cj = np.nonzero(coast)
    cl = np.cos(np.radians(LAT[ci]))
    R = 4
    src, dst, wkm = [], [], []
    for di in range(-R, R + 1):
        for dj in range(-R - 2, R + 3):
            if max(abs(di), abs(dj)) <= 1: continue
            i2, j2 = ci + di, cj + dj
            ok = (i2 >= 0) & (i2 < NY) & (j2 >= 0) & (j2 < NX)
            i2c, j2c = np.clip(i2, 0, NY - 1), np.clip(j2, 0, NX - 1)
            ok &= coast[i2c, j2c]
            dist = CELL_KM * np.hypot(di, dj * cl)
            ok &= dist <= max_km
            ns = int(np.ceil(CELL_KM * np.hypot(di, dj) / 4)) + 1
            fr = np.linspace(0, 1, ns)
            ii = np.rint(ci[:, None] + di * fr[None]).astype(int)
            jj = np.rint(cj[:, None] + dj * fr[None]).astype(int)
            ii = np.clip(ii, 0, NY - 1); jj = np.clip(jj, 0, NX - 1)
            wf = water[ii, jj].mean(axis=1)
            ok &= wf >= 0.25
            src.append((ci * NX + cj)[ok]); dst.append((i2c * NX + j2c)[ok]); wkm.append(np.maximum(dist * wf - 25, 2)[ok])
    return np.concatenate(src), np.concatenate(dst), np.concatenate(wkm).astype(np.float32)

def build_crossings():
    ptr = np.zeros((len(SL_BINS), N + 1), np.int64); idx = []; w = []; base = 0
    for b, sl in enumerate(SL_BINS):
        s, d, k = _crossing_edges(water_mask(sl))
        o = np.argsort(s, kind='stable'); s, d, k = s[o], d[o], k[o]
        counts = np.bincount(s, minlength=N)
        ptr[b, 1:] = base + np.cumsum(counts); ptr[b, 0] = base
        idx.append(d); w.append(k); base += len(d)
    return ptr, np.concatenate(idx).astype(np.int64), np.concatenate(w).astype(np.float32)

# ---------------------------------------------------------------- climate indices
def climate(bp, p):
    """Returns (H: Saharo-Arabian humidity 0..1, M: monsoon index 0..1) at years BP."""
    ka = bp / 1000
    H = np.full_like(ka, p['H_bg'], dtype=float)
    def pulse(a, b, amp):           # smooth-edged window a>b in ka
        return amp * (1 / (1 + np.exp(-(a - ka) / 0.6))) * (1 / (1 + np.exp(-(ka - b) / 0.6)))
    H = np.maximum(H, pulse(p['h5c_a'], p['h5c_b'], p['h5c_amp']))
    H = np.maximum(H, pulse(p['h5a_a'], p['h5a_b'], p['h5a_amp']))
    H = np.maximum(H, pulse(p['h3_a'], p['h3_b'], p['h3_amp']))
    mis4 = (ka < p['h5a_b']) & (ka > 60)
    H = np.where(mis4, np.minimum(H, p['H_mis4']), H)
    M = np.where(ka > p['h5a_b'], 0.95, np.where(ka > 60, p['M_mis4'], p['M_mis3']))
    return H, M

# ---------------------------------------------------------------- strategy speed surfaces
def speed_surfaces(p, strategy, noise):
    """speed[K, N] in km/yr; 0 = uninhabitable at that time."""
    sp = np.zeros((K, N), np.float32); coastk = np.zeros((K, N), np.bool_)
    H, M = climate(T_K.astype(float), p)
    for k in range(K):
        sl = sea_level(T_K[k], p['sl_off'], p['sl_jit'])
        water = water_mask(sl)
        land = ~water
        coast = land & ndimage.binary_dilation(water & ~LAKE, structure=np.ones((3, 3)))
        h, m = H[k], M[k]
        glacial = np.clip(-sl / 100, 0, 1)
        # --- terrestrial (big-game) habitat quality
        sT = np.select([BIOME == 0, BIOME == 1, BIOME == 2, BIOME == 3, BIOME == 4],
                       [0.02 + 0.85 * h, 0.3 + 0.45 * h, 0.45 + 0.4 * m, 0.25 + 0.35 * glacial, 0.8])
        sT = ndimage.gaussian_filter(sT, 1.5)          # soften hand-drawn biome edges
        sT = np.where(RIVER, np.maximum(sT, 0.6), sT)
        sT = np.where(GULF_BASIN, np.maximum(sT, 0.45), sT)       # rivers + springs on the exposed Gulf floor
        sT = sT * ALT_PEN
        # --- coastal habitat quality (shellfish, reef fish, mangrove, strandloping)
        fw = np.where(BIOME == 0, p['fw_floor'] + (1 - p['fw_floor']) * h,
             np.where(BIOME == 1, 0.5 * (1 + p['fw_floor']) + 0.5 * (1 - p['fw_floor']) * h, 1.0))
        fw = ndimage.gaussian_filter(fw, 1.5)
        sC = p['c_adv'] * fw
        if strategy == 'coastal':
            s = np.where(coast, sC, p['inland_frac'] * sT)
        elif strategy == 'terrestrial':
            s = sT
        else:  # mixed generalist
            s = np.where(coast, np.maximum(sC, sT), sT)
        s = s * noise
        s = np.where(land & (s >= 0.12), s, 0.0)
        sp[k] = (p['v_max'] * s).ravel(); coastk[k] = coast.ravel()
    return sp, coastk

# ---------------------------------------------------------------- time-dependent Dijkstra (numba)
@nb.njit(cache=True)
def _dijkstra(speed, coastk, biome, scode, pen, nbr_d, offs, ny, nx, seeds, tau0, slbin, cptr, cidx, cw, lam, w0, dmax, dt, tau_max):
    Kk, Nn = speed.shape
    arr = np.full(Nn, np.inf); pred = np.full(Nn, -1, np.int64); via_sea = np.zeros(Nn, np.bool_)
    heap = [(0.0, np.int64(0))]; heap.pop()
    for s in seeds:
        k0 = int(tau0 // dt)
        if speed[k0, s] > 0:
            arr[s] = tau0; heap.append((tau0, s))
    heapq.heapify(heap)
    while len(heap) > 0:
        t, u = heapq.heappop(heap)
        if t > arr[u] or t > tau_max: continue
        k0 = int(t // dt)
        if k0 >= Kk: continue
        r = u // nx; c = u % nx
        for e in range(8):
            r2 = r + offs[e, 0]; c2 = c + offs[e, 1]
            if r2 < 0 or r2 >= ny or c2 < 0 or c2 >= nx: continue
            v = r2 * nx + c2
            k = k0
            while k < Kk:
                su = speed[k, u]
                if su <= 0: break                      # population at u is gone
                sv = speed[k, v]
                if sv > 0:
                    dep = max(t, k * dt)
                    na = dep + nbr_d[r, e] / (0.5 * (su + sv))
                    # learning cost when entering a new biome, unless a coast-user stays on the shore
                    if biome[u] != biome[v] and not (scode != 1 and coastk[k, u] and coastk[k, v]):
                        na += pen
                    if na < arr[v]:
                        arr[v] = na; pred[v] = u; via_sea[v] = False
                        heapq.heappush(heap, (na, v))
                    break
                k += 1
        b = slbin[k0]
        if speed[k0, u] > 0:
            for q in range(cptr[b, u], cptr[b, u + 1]):
                v = cidx[q]; w = cw[q]
                if w > dmax or speed[k0, v] <= 0: continue
                na = t + w0 * np.exp(w / lam)
                if na < arr[v]:
                    arr[v] = na; pred[v] = u; via_sea[v] = True
                    heapq.heappush(heap, (na, v))
    return arr, pred, via_sea

def run(p, strategy, seeds, t0_bp, noise, X, tau_max_bp=T_END):
    sp, coastk = speed_surfaces(p, strategy, noise)
    scode = {'coastal': 0, 'terrestrial': 1, 'mixed': 2}[strategy]
    slb = np.array([np.argmin(np.abs(SL_BINS - sea_level(T_K[k], p['sl_off'], p['sl_jit']))) for k in range(K)])
    arr, pred, via = _dijkstra(sp, coastk, BIOME.ravel(), scode, float(p['biome_pen']), NBR_D, OFFS, NY, NX, np.asarray(seeds, np.int64), float(bp_to_tau(t0_bp)),
                               slb, X[0], X[1], X[2], p['lam'], p['w0'], p['dmax'], float(DT),
                               float(bp_to_tau(tau_max_bp)))
    return tau_to_bp(arr), pred, via

# ---------------------------------------------------------------- parameter sampler
def sample_params(rng):
    lu = lambda a, b: float(np.exp(rng.uniform(np.log(a), np.log(b))))
    return dict(
        sl_off=rng.normal(0, 8), sl_jit=rng.uniform(-1500, 1500),
        H_bg=rng.uniform(0.0, 0.15),
        h5c_a=rng.uniform(102, 105), h5c_b=rng.uniform(95, 98), h5c_amp=rng.uniform(0.5, 1.0),
        h5a_a=rng.uniform(82, 86), h5a_b=rng.uniform(70, 73), h5a_amp=rng.uniform(0.5, 1.0),
        h3_a=rng.uniform(57, 60), h3_b=rng.uniform(52, 55), h3_amp=rng.uniform(0.1, 0.5),
        H_mis4=rng.uniform(0.0, 0.12),
        M_mis4=rng.uniform(0.2, 0.6), M_mis3=rng.uniform(0.5, 0.85),
        v_max=lu(0.7, 6.0),                 # km/yr front speed in ideal habitat
        c_adv=rng.uniform(0.7, 1.2),        # coastal habitat quality relative to best savanna
        fw_floor=rng.uniform(0.05, 0.7),    # freshwater on hyper-arid coasts in dry phases (springs)
        inland_frac=rng.uniform(0.1, 0.35), # how well coast-specialists cope away from the coast
        biome_pen=rng.uniform(0, 1500),     # years to re-tool when entering a new biome
        lam=rng.uniform(15, 45), w0=lu(20, 500), dmax=rng.uniform(60, 130),
    )

def noise_field(rng, sd=0.25, smooth=3):
    f = ndimage.gaussian_filter(rng.normal(size=(NY, NX)), smooth)
    f = f / f.std() * sd
    return np.exp(f - sd ** 2 / 2)

def cell(lat, lon):
    i = int(np.argmin(np.abs(LAT - lat))); j = int(np.argmin(np.abs(LON - lon)))
    return i * NX + j
