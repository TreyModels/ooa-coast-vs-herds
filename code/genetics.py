"""
Genetic consequences of each dispersal strategy.

Wave 2 = the ancestors of everyone living outside Africa: spreads from the Persian-Plateau hub
after the shared Neanderthal admixture (start 50-46 ka). Its colonisation tree (who founded whom)
gives (a) population split times, (b) serial-founder loss of diversity, scaled by how narrow the
expanding front was, and (c) how much ancestry it picks up from an EARLIER wave (wave 1 = the
long-chronology runs, which may already be sitting along the route) as it moves through their land.
"""
import os
os.environ['OOA_T_END'] = '25000'
import numpy as np, model as M, json
from multiprocessing import Pool
C = np.load('cross.npz'); X = (C['ptr'], C['idx'], C['w'])
STRATS = ['coastal', 'terrestrial', 'mixed']
HUB = np.flatnonzero((M._box(46, 58, 25.5, 35) & ~M.water_mask(-60) & (M.ZMED < 2200)).ravel())
POPS = {  # present-day / ancient reference populations
    'Levant (West Eurasian side)': (32.5, 35.5),
    'Indus (Pakistan)': (27.0, 68.5),
    'South India (AASI)': (11.0, 77.0),
    'Andamans (Onge)': (12.0, 92.8),
    'Malay Peninsula (Semang)': (5.0, 101.5),
    'North China (Tianyuan)': (39.7, 115.9),
    'Luzon (Ayta)': (15.2, 120.5),
    'New Guinea (Papuans)': (-5.8, 144.5),
    'N Australia (Aboriginal)': (-13.0, 133.5),
}
def near(lat, lon, arr):
    i = int(np.argmin(abs(M.LAT - lat))); j = int(np.argmin(abs(M.LON - lon)))
    best = None
    for r in range(0, 5):
        sl = (slice(max(0, i - r), i + r + 1), slice(max(0, j - r), j + r + 1))
        idx = np.ravel_multi_index(np.mgrid[sl].reshape(2, -1), (M.NY, M.NX))
        ok = idx[np.isfinite(arr[idx])]
        if len(ok): return int(ok[np.argmax(arr[ok])])
    return None

def front_width(arr, cells, win=400, R=2):
    A = arr.reshape(M.NY, M.NX); out = []
    for c in cells:
        i, j = divmod(int(c), M.NX)
        blk = A[max(0, i - R):i + R + 1, max(0, j - R):j + R + 1]
        out.append(int(np.sum(np.isfinite(blk) & (np.abs(blk - A[i, j]) <= win))))
    return np.array(out)

def one(i):
    rng = np.random.default_rng(50_000 + i)
    p = M.sample_params(rng); noise = M.noise_field(rng)
    t0 = rng.uniform(46_000, 50_000)
    res = {'t0': t0, 'p': p}
    for s in STRATS:
        arr, pred, via = M.run(p, s, HUB, t0, noise, X, tau_max_bp=25_000)
        d = {}
        for name, (la, lo) in POPS.items():
            c = near(la, lo, arr)
            if c is None: d[name] = None; continue
            path = []; u = c
            while u >= 0: path.append(u); u = pred[u]
            path = np.array(path[::-1])                     # hub -> population
            d[name] = dict(path=path, t=arr[path].astype(np.float32), w=front_width(arr, path))
        res[s] = d
    return res

if __name__ == '__main__':
    n = 150
    with Pool(os.cpu_count()) as pool: R = pool.map(one, range(n), chunksize=5)
    import pickle; pickle.dump(R, open('wave2.pkl', 'wb'))
    print('done', sum(r['mixed']['New Guinea (Papuans)'] is not None for r in R), 'of', n, 'reach PNG (mixed)')
