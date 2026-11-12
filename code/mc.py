"""Paired Monte Carlo: same parameter draw + same noise field for every strategy."""
import numpy as np, sys, os, time, model as M
from multiprocessing import Pool

C = np.load('cross.npz'); X = (C['ptr'], C['idx'], C['w'])
STRATS = ['coastal', 'terrestrial', 'mixed']
PRESENT_WATER = M.water_mask(0.0).ravel()
SAHUL = (((M.LON2 >= 131.5) & (M.LAT2 <= 0)) | ((M.LAT2 <= -11.2) & (M.LON2 >= 121))).ravel()
HORN = np.flatnonzero((M._box(36, 43.5, 8, 15.5) & ~M.water_mask(-40) & (M.ZMED < 2000)).ravel())
NEAREAST = np.flatnonzero((M._box(34.5, 48, 29.5, 34) & ~M.water_mask(-40) & (M.ZMED < 1500)).ravel())

def region(lat, lon, r=2):
    i = int(np.argmin(abs(M.LAT - lat))); j = int(np.argmin(abs(M.LON - lon)))
    m = np.zeros((M.NY, M.NX), bool); m[max(0, i - r):i + r + 1, max(0, j - r):j + r + 1] = True
    return np.flatnonzero(m.ravel())
SITES = {
    'S Arabia coast (Yemen)': region(13.0, 45.0),
    'Jebel Faya (UAE)': region(25.1, 55.8),
    'Nefud palaeolakes (N Arabia)': region(28.0, 40.5, 3),
    'Levant (Skhul/Manot)': region(32.6, 35.4),
    'Jwalapuram (S India)': region(15.3, 78.1),
    'Lida Ajer (Sumatra)': region(-0.25, 100.3),
    'Sulawesi (Maros / Muna)': np.r_[region(-5.0, 119.7), region(-4.9, 122.6)],
    'Madjedbebe (N Australia)': region(-12.5, 132.9),
}

def route(pred, end):
    path = [];
    u = end
    while u >= 0 and len(path) < 100000:
        path.append(u); u = pred[u]
    return np.array(path)

def one(args):
    scen, i = args
    rng = np.random.default_rng(10_000 * (scen == 'S') + i)
    p = M.sample_params(rng)
    noise = M.noise_field(rng)
    if scen == 'L':
        t0 = rng.uniform(72_000, 85_000); seeds = HORN
    else:
        t0 = rng.uniform(47_000, 51_000); seeds = NEAREAST
    out = dict(scen=scen, i=i, t0=t0, **p)
    routes = {}
    for s in STRATS:
        arr, pred, via = M.run(p, s, seeds, t0, noise, X)
        sah = np.where(SAHUL & np.isfinite(arr), arr, -np.inf)
        first = int(np.argmax(sah)); t_sahul = sah[first]
        out[f'{s}|Sahul'] = t_sahul
        for name, idx in SITES.items():
            out[f'{s}|{name}'] = np.max(np.where(np.isfinite(arr[idx]), arr[idx], -np.inf))
        if np.isfinite(t_sahul):
            path = route(pred, first)
            out[f'{s}|route_cells'] = len(path)
            out[f'{s}|route_sunk'] = PRESENT_WATER[path].mean()
            out[f'{s}|route_crossings'] = int(via[path].sum())
            routes[s] = path
        else:
            out[f'{s}|route_sunk'] = np.nan
        # occupied-by-snapshot maps (for probability maps)
        snap = 65_000 if scen == 'L' else 45_000
        routes[s + '_occ'] = np.flatnonzero(arr >= snap)
    return out, routes

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    jobs = [(sc, i) for sc in ['L', 'S'] for i in range(n)]
    t = time.time()
    with Pool(os.cpu_count()) as pool:
        res = pool.map(one, jobs, chunksize=4)
    print('done', time.time() - t)
    import pandas as pd
    df = pd.DataFrame([r[0] for r in res]); df.to_pickle('mc.pkl')
    dens = {}
    for (sc, i), (_, routes) in zip(jobs, res):
        for k, v in routes.items():
            key = f'{sc}|{k}'
            if key not in dens: dens[key] = np.zeros(M.N, np.int32)
            dens[key][v] += 1
    np.savez_compressed('density.npz', n=n, **dens)
