"""Paired Monte Carlo, two climates: hand-drawn boxes ('box', = original mc.py) vs LOVECLIM ('ibs').

For run (scen, i) the base parameters, noise field and start time are drawn exactly as in mc.py, so the
'box' rows reproduce the original runs. LOVECLIM-only parameters come from a separate RNG.
Every strategy x climate combination uses the same draw.

Run from a work folder next to code/ (needs cross.npz from build.py):
    python ../code/mc_ibs.py 400          # n runs per chronology
    python ../code/mc_ibs.py 4 L          # quick smoke test, long chronology only
Writes mc_ibs.pkl and density_ibs.npz.
"""
import numpy as np, sys, os, time
from multiprocessing import Pool
import model as M, model_ibs as MI, climate_ibs as CI
import mc as MC                                 # SITES, SAHUL, HORN, NEAREAST, X, route()

STRATS = MC.STRATS
CLIMS = {'box': M.run, 'ibs': MI.run}
if os.environ.get('OOA_CLIMS'): CLIMS = {k: CLIMS[k] for k in os.environ['OOA_CLIMS'].split(',')}
C_MIX = float(os.environ.get('OOA_C_MIX', 0.0))       # review: cost of the mixed strategy
TAG = os.environ.get('OOA_TAG', '')


def one(args):
    scen, i = args
    rng = np.random.default_rng(10_000 * (scen == 'S') + i)        # identical to mc.py
    p = M.sample_params(rng)
    noise = M.noise_field(rng)
    if scen == 'L':
        t0 = rng.uniform(72_000, 85_000); seeds = MC.HORN
    else:
        t0 = rng.uniform(47_000, 51_000); seeds = MC.NEAREAST
    p.update(MI.sample_extra(np.random.default_rng(900_000 + 10_000 * (scen == 'S') + i)))
    p['c_mix'] = C_MIX
    out = dict(scen=scen, i=i, t0=t0, **p)
    routes = {}
    snap = 65_000 if scen == 'L' else 45_000
    for cl, runner in CLIMS.items():
        for s in STRATS:
            arr, pred, via = runner(p, s, seeds, t0, noise, MC.X)
            sah = np.where(MC.SAHUL & np.isfinite(arr), arr, -np.inf)
            first = int(np.argmax(sah)); t_sahul = sah[first]
            pre = f'{cl}|{s}'
            out[f'{pre}|Sahul'] = t_sahul
            for name, idx in MC.SITES.items():
                out[f'{pre}|{name}'] = np.max(np.where(np.isfinite(arr[idx]), arr[idx], -np.inf))
            if np.isfinite(t_sahul):
                path = MC.route(pred, first)
                out[f'{pre}|route_cells'] = len(path)
                out[f'{pre}|route_sunk'] = MC.PRESENT_WATER[path].mean()
                out[f'{pre}|route_crossings'] = int(via[path].sum())
                routes[pre] = path
            else:
                out[f'{pre}|route_sunk'] = np.nan
            routes[pre + '_occ'] = np.flatnonzero(arr >= snap)
    return out, routes


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    scens = list(sys.argv[2]) if len(sys.argv) > 2 else ['L', 'S']
    CI.load()                                          # build the climate cache once before forking
    jobs = [(sc, i) for sc in scens for i in range(n)]
    t = time.time()
    with Pool(os.cpu_count()) as pool:
        res = []
        for k, r in enumerate(pool.imap(one, jobs, chunksize=2), 1):
            res.append(r)
            if k % max(1, len(jobs) // 20) == 0 or k == len(jobs):
                print(f'{k}/{len(jobs)} runs  {time.time() - t:.0f}s', flush=True)
    import pandas as pd
    df = pd.DataFrame([r[0] for r in res]); df.to_pickle(f'mc_ibs{TAG}.pkl')
    dens = {}
    for (sc, i), (_, routes) in zip(jobs, res):
        for k, v in routes.items():
            key = f'{sc}|{k}'
            if key not in dens: dens[key] = np.zeros(M.N, np.int32)
            dens[key][v] += 1
    np.savez_compressed(f'density_ibs{TAG}.npz', n=n, **dens)
    print('done', round(time.time() - t), 's;', len(df), 'runs ->', os.path.abspath(f'mc_ibs{TAG}.pkl'))
