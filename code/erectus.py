"""
Homo erectus-grade dispersal, 2.6 Ma -> 0.1 Ma, on the same paleo-landscape engine.

3 foraging strategies x 3 sea-crossing abilities:
  walk  - land bridges only (gaps up to ~5 km after grid correction)
  raft  - rare accidental crossings (natural rafts, tsunami debris): tens to hundreds of kyr per crossing
  canoe - deliberate dugout-canoe crossings, same ability range as the sapiens runs
Two sea-level histories, sampled 50/50:
  classic - shallower 41-kyr lowstands before the Mid-Pleistocene Transition
  clark25 - Clark et al. 2025: early lowstands comparable to the Last Glacial Maximum
Scored on WHERE hominins are / are not found, with generous date limits.
"""
import os
import numpy as np, model as M
from scipy import ndimage
from multiprocessing import Pool

C = np.load('cross.npz'); X = (C['ptr'], C['idx'], C['w'])
T0_MAX, T_END, DT = 2_700_000, 100_000, 5000
K = (T0_MAX - T_END) // DT + 1
T_K = T0_MAX - np.arange(K) * DT
STRATS = ['coastal', 'terrestrial', 'mixed']
CROSS = ['walk', 'raft', 'canoe']
SEEDS = np.flatnonzero((M._box(35, 43.5, 2, 15) & ~M.water_mask(-40) & (M.ZMED < 2200)).ravel())

def region(lat, lon, r=2):
    i = int(np.argmin(abs(M.LAT - lat))); j = int(np.argmin(abs(M.LON - lon)))
    m = np.zeros((M.NY, M.NX), bool); m[max(0, i - r):i + r + 1, max(0, j - r):j + r + 1] = True
    return np.flatnonzero(m.ravel())
# (lat, lon, minimum age in years, radius) – found
PRESENT = {
    'Shangchen, Loess Plateau (2.1 Ma)': (34.1, 109.3, 2_120_000, 2),
    'Dmanisi, Georgia (1.85 Ma)': (41.3, 44.3, 1_850_000, 2),
    'Yuanmou, S China (1.7 Ma)': (25.7, 101.9, 1_700_000, 2),
    'Nihewan, N China (1.66 Ma)': (40.2, 114.6, 1_660_000, 2),
    'Sangiran, Java (1.5 Ma)': (-7.45, 110.8, 1_500_000, 2),
    'Attirampakkam, India (1.5 Ma)': (13.2, 79.9, 1_500_000, 2),
    'Ubeidiya, Levant (1.4 Ma)': (32.7, 35.55, 1_400_000, 2),
    'Calio, Sulawesi (1.04 Ma)': (-4.4, 119.9, 1_040_000, 2),
    'Wolo Sege, Flores (1.02 Ma)': (-8.7, 121.1, 1_020_000, 2),
    'Kalinga, Luzon (0.71 Ma)': (17.6, 121.4, 709_000, 2),
}
# never found before modern humans
ABSENT = {
    'Sahul (New Guinea + Australia)': np.flatnonzero((((M.LON2 >= 131.5) & (M.LAT2 <= 0)) | ((M.LAT2 <= -11.2) & (M.LON2 >= 121))).ravel()),
    'Andaman Islands': region(12.0, 92.8, 2),
    'Socotra': region(12.5, 53.9, 2),
    'Japan (Honshu)': region(36.0, 138.0, 4),
}
PRES_IDX = {k: region(v[0], v[1], v[3]) for k, v in PRESENT.items()}

def climate_path(p):
    t = T_K.astype(float)
    g41 = 0.5 * (1 + np.sin(2 * np.pi * t / 41_000 + p['phi']))
    saw = ((p['phi2'] * 1e4 - t) % 100_000) / 100_000            # slow build-up, fast termination
    g100 = np.clip(saw ** 1.5 * 1.1, 0, 1)
    w = np.clip((1_200_000 - t) / 300_000, 0, 1)                   # blend across the MPT
    g = (1 - w) * g41 + w * np.maximum(g100, 0.6 * g41)
    if p['sl_model'] == 'clark25':
        low = np.full_like(t, p['low_late'])
    else:
        low = np.where(t > 1_200_000, p['low_early'], np.where(t < 900_000, p['low_late'],
              p['low_early'] + (p['low_late'] - p['low_early']) * (1_200_000 - t) / 300_000))
    sl = -low * g ** 1.3
    H = p['h_amp'] * (1 - g) ** 1.5                                 # Arabia wet in interglacials
    Mo = 0.45 + 0.5 * (1 - g)
    return g, sl, H, Mo

# Optional: real paleoclimate habitat (Timmermann et al. 2022) instead of the hand-drawn biome boxes.
HAB = None
if os.environ.get('ERECTUS_HABITAT'):
    from climate_inputs import HabitatField
    HAB = HabitatField(os.environ['ERECTUS_HABITAT'], os.environ.get('ERECTUS_HABITAT_VAR', 'erectus'))

_MASK = {}
def masks(b):
    if b not in _MASK:
        water = M.water_mask(float(M.SL_BINS[b])); land = ~water
        coast = land & ndimage.binary_dilation(water & ~M.LAKE, structure=np.ones((3, 3)))
        _MASK[b] = (water, land, coast)
    return _MASK[b]

def surfaces(p, strategy, noise):
    g, sl, H, Mo = climate_path(p)
    slb = np.array([int(np.argmin(np.abs(M.SL_BINS - max(s, -120)))) for s in sl])
    hq = np.round(H * 5) / 5; mq = np.round(Mo * 5) / 5; gq = np.round(g * 4) / 4
    sp = np.zeros((K, M.N), np.float32); ck = np.zeros((K, M.N), np.bool_); cache = {}
    cold = np.clip((M.LAT2 - 36) / 10, 0, 1)
    for k in range(K):
        key = (slb[k], hq[k], mq[k], gq[k]) if HAB is None else (slb[k], hq[k], mq[k], gq[k], int(T_K[k] // 25_000))
        if key not in cache:
            water, land, coast = masks(slb[k]); h, m, gg = key[1], key[2], key[3]
            B = M.BIOME
            sT = np.select([B == 0, B == 1, B == 2, B == 3, B == 4],
                           [0.02 + 0.85 * h, 0.3 + 0.45 * h, 0.45 + 0.4 * m, 0.3 + 0.3 * gg, 0.8])
            sT = ndimage.gaussian_filter(sT, 1.5)
            sT = np.where(M.RIVER, np.maximum(sT, 0.6), sT); sT = np.where(M.GULF_BASIN, np.maximum(sT, 0.45), sT)
            sT = sT * M.ALT_PEN * (1 - p['cold'] * gg * cold)       # glacial cold limits the north
            if HAB is not None:                                      # replace with simulated habitat suitability
                sT = HAB.at(float(T_K[k])) * M.ALT_PEN
                h = float(np.mean(sT[M.BIOME == 0])) if (M.BIOME == 0).any() else h
            fw = np.where(B == 0, p['fw_floor'] + (1 - p['fw_floor']) * h,
                 np.where(B == 1, 0.5 * (1 + p['fw_floor']) + 0.5 * (1 - p['fw_floor']) * h, 1.0))
            fw = ndimage.gaussian_filter(fw, 1.5)
            sC = p['c_adv'] * fw * (1 - p['cold'] * gg * cold)
            if strategy == 'coastal': s = np.where(coast, sC, p['inland_frac'] * sT)
            elif strategy == 'terrestrial': s = sT
            else: s = np.where(coast, np.maximum(sC, sT), sT)
            s = s * noise
            s = np.where(land & (s >= 0.12), s, 0.0)
            cache[key] = ((p['v_max'] * s).ravel().astype(np.float32), coast.ravel())
        sp[k], ck[k] = cache[key]
    return sp, ck, slb

def sample(rng):
    lu = lambda a, b: float(np.exp(rng.uniform(np.log(a), np.log(b))))
    return dict(phi=rng.uniform(0, 2 * np.pi), phi2=rng.uniform(0, 10),
                sl_model='clark25' if rng.random() < 0.5 else 'classic',
                low_early=rng.uniform(50, 85), low_late=rng.uniform(100, 125),
                h_amp=rng.uniform(0.4, 1.0), cold=rng.uniform(0.3, 0.9),
                v_max=lu(0.1, 2.0), c_adv=rng.uniform(0.7, 1.2), fw_floor=rng.uniform(0.05, 0.7),
                inland_frac=rng.uniform(0.1, 0.35), biome_pen=rng.uniform(0, 5000),
                t0=rng.uniform(2_250_000, 2_600_000),
                raft_w0=lu(20_000, 300_000), raft_lam=rng.uniform(15, 30), raft_dmax=rng.uniform(30, 80),
                canoe_w0=lu(50, 1000), canoe_lam=rng.uniform(15, 40), canoe_dmax=rng.uniform(40, 100))

def run_one(i):
    rng = np.random.default_rng(900_000 + i); p = sample(rng); noise = M.noise_field(rng)
    out = {'i': i, **{k: v for k, v in p.items()}}; reach = {}
    for s in STRATS:
        sp, ck, slb = surfaces(p, s, noise)
        scode = {'coastal': 0, 'terrestrial': 1, 'mixed': 2}[s]
        for c in CROSS:
            if c == 'walk': w0, lam, dmax = 1.0, 1.0, 5.0
            else: w0, lam, dmax = p[c + '_w0'], p[c + '_lam'], p[c + '_dmax']
            arr, pred, via = M._dijkstra(sp, ck, M.BIOME.ravel(), scode, float(p['biome_pen']), M.NBR_D, M.OFFS,
                                         M.NY, M.NX, SEEDS, float(T0_MAX - p['t0']), slb, X[0], X[1], X[2],
                                         lam, w0, dmax, float(DT), float(T0_MAX - T_END))
            bp = T0_MAX - arr            # years BP; -inf if never
            bp = np.where(bp >= T_END, bp, -np.inf)
            for name, (la, lo, age, r) in PRESENT.items():
                v = bp[PRES_IDX[name]]; out[f'{s}|{c}|{name}'] = float(np.max(np.where(np.isfinite(v), v, -np.inf)))
            for name, idx in ABSENT.items():
                v = bp[idx]; out[f'{s}|{c}|NOT {name}'] = float(np.max(np.where(np.isfinite(v), v, -np.inf)))
            reach[f'{s}|{c}'] = np.isfinite(bp) & (bp >= T_END)
    return out, reach

if __name__ == '__main__':
    import sys, pandas as pd
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(os.cpu_count()) as pool: R = pool.map(run_one, range(n), chunksize=2)
    pd.DataFrame([r[0] for r in R]).to_pickle('erectus.pkl')
    dens = {k: np.zeros(M.N, np.int32) for k in R[0][1]}
    for _, rch in R:
        for k, v in rch.items(): dens[k] += v
    np.savez_compressed('erectus_reach.npz', n=n, **dens)
    print('done', n)
