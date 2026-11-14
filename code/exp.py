"""Experiments 1 (rejection-ABC) and 2 (drowned-site hotspots) on one mc run.

usage: python exp.py TAG [--no-nefud]
  reads  mc_TAG.pkl, density_TAG.npz           (written by mc.py --tag TAG)
  writes exp_TAG[_nonefud].json, map_drowned_TAG[_nonefud].png
--no-nefud drops the Nefud >=70 ka term from the ABC acceptance (review item 3).
"""
import numpy as np, pandas as pd, json, sys, model as M, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb

TAG = sys.argv[1] if len(sys.argv) > 1 else 'old'
NO_NEFUD = '--no-nefud' in sys.argv
SUF = f"{TAG}_nonefud" if NO_NEFUD else TAG

df = pd.read_pickle(f'mc_{TAG}.pkl'); L = df[df.scen == 'L']; S = ['coastal', 'terrestrial', 'mixed']
out = {'tag': TAG, 'no_nefud': NO_NEFUD}
# ---- Experiment 1: rejection-ABC. Keep draws whose outcome matches the evidence.
ev = {}
for s in S:
    ok = (L[f'{s}|Sulawesi (Maros / Muna)'] >= 67800) & (L[f'{s}|Lida Ajer (Sumatra)'] >= 63000)
    if not NO_NEFUD:
        ok &= (L[f'{s}|Nefud palaeolakes (N Arabia)'] >= 70000)
    ev[s] = ok
tot = sum(ev[s].sum() for s in S)
out['n_accepted'] = {s: int(ev[s].sum()) for s in S}
out['post_strategy'] = {s: round(float(ev[s].sum() / tot), 3) for s in S}
def qs(x): return [round(float(np.percentile(x, q)), 2) for q in (10, 50, 90)]
acc = pd.concat([L[ev[s]].assign(strat=s) for s in S])
out['abc'] = {}
for prm, scale in [('v_max', 1), ('t0', 1 / 1000), ('c_adv', 1), ('fw_floor', 1), ('biome_pen', 1), ('lam', 1), ('w0', 1)]:
    out['abc'][prm] = {'prior': qs(L[prm] * scale), 'post_all': qs(acc[prm] * scale), **{s: qs(L[ev[s]][prm] * scale) for s in S}}
# coastal-vs-terrestrial head-to-head among accepted mixed-evidence draws (same draws)
ok = ev['mixed']; c = L['coastal|Sahul'].replace(-np.inf, 0); t = L['terrestrial|Sahul'].replace(-np.inf, 0)
out['h2h_given_evidence'] = round(float((c > t)[ok].mean()), 3)
# speed histogram prior vs posterior (mixed)
bins = np.exp(np.linspace(np.log(0.7), np.log(6), 9))
out['vhist'] = {'edges': [round(b, 2) for b in bins], 'prior': np.histogram(L.v_max, bins)[0].tolist(), 'post': np.histogram(acc.v_max, bins)[0].tolist(),
  'rate': [round(float(ok[(L.v_max >= a) & (L.v_max < b)].mean()), 3) for a, b in zip(bins[:-1], bins[1:])]}
out['t0rate'] = [{'bin': f'{a}–{a+2.5:g}', 'rate': round(float(ok[(L.t0 >= a * 1000) & (L.t0 < (a + 2.5) * 1000)].mean()), 3)} for a in [72, 74.5, 77, 79.5, 82]]

# ---- Experiment 2: where to look for drowned sites
D = np.load(f'density_{TAG}.npz'); n = int(D['n'])
dens = (D['L|coastal'] + D['L|mixed']) / (2 * n)
seabed = M.water_mask(0.0).ravel(); depth = -M.ZMED.ravel()
mask = seabed & (dens >= 0.05) & (depth > 0) & (depth <= 130)
# shelf area above -120 m from the model grid: 20' cells (CELL_KM square, lon width scaled by cos lat)
cell_km2 = (M.CELL_KM ** 2) * np.cos(np.radians(M.LAT2)).ravel()
shelf = seabed & (depth > 0) & (depth <= 120)
regions = {'Red Sea & Bab-el-Mandeb': (32, 44, 11, 30), 'Gulf of Aden & S Arabia': (43, 53, 10, 17), 'Oman & Persian Gulf': (48, 60, 17, 31),
 'Makran & Indus delta': (60, 70, 22, 27), 'W India shelf (Gujarat–Kerala)': (68, 77.5, 7, 24), 'Sri Lanka & Palk Strait': (77.5, 82.5, 5, 11),
 'Bay of Bengal & Myanmar': (80, 98, 8, 23), 'Sunda shelf (Malaya, Sumatra, Java)': (95, 116, -9, 8), 'Wallacea': (116, 131, -11, 2), 'Sahul shelf': (121, 142, -15, -8)}
rows = []; LON = M.LON2.ravel(); LAT = M.LAT2.ravel(); W = dens * mask
for name, (a, b, c_, d) in regions.items():
    box = (LON >= a) & (LON <= b) & (LAT >= c_) & (LAT <= d)
    m = box & mask
    shelf_k = float(cell_km2[box & shelf].sum())
    if m.sum() == 0: continue
    score = float(W[m].sum())
    rows.append(dict(region=name, score=round(score, 3), share=round(score / float(W.sum()), 3), depth_med=int(np.median(depth[m])),
       shallow=round(float(W[m & (depth <= 40)].sum() / W[m].sum()), 2),
       shelf_km2_above_120m=round(shelf_k), score_per_1000km2_shelf=round(score / (shelf_k / 1000), 5) if shelf_k > 0 else None))
rows.sort(key=lambda r: -r['share']); out['hotspots'] = rows
rows_n = sorted(rows, key=lambda r: -(r['score_per_1000km2_shelf'] or 0)); out['hotspots_normalised_rank'] = [r['region'] for r in rows_n]
out['drowned_depth'] = [round(float(np.percentile(depth[mask], q, method='inverted_cdf')), 0) for q in (10, 25, 50, 75, 90)]
json.dump(out, open(f'exp_{SUF}.json', 'w'), indent=1); print(json.dumps(out, indent=1))
# map
r = (M.LAT <= 32) & (M.LAT >= -15); c = (M.LON >= 32) & (M.LON <= 142); ext = [M.LON[c][0], M.LON[c][-1], M.LAT[r][-1], M.LAT[r][0]]
land_now = ~M.water_mask(0.0); land_lo = ~M.water_mask(-85.0)
rgb = np.ones(land_now.shape + (3,)) * np.array(to_rgb('#e9eef3')); rgb[land_lo] = to_rgb('#d8d3c8'); rgb[land_now] = to_rgb('#ece8e0')
fig = plt.figure(figsize=(9, 3.9)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
ax.imshow(rgb[np.ix_(r, c)], extent=ext, interpolation='nearest')
dm = np.where(mask, depth, np.nan).reshape(M.NY, M.NX)
cm = LinearSegmentedColormap.from_list('d', ['#9ec5f0', '#2a78d6', '#0d3a75'])
ii = np.flatnonzero(mask)
ax.scatter(LON[ii], LAT[ii], s=6 + 40 * np.clip(dens[ii] / 0.5, 0, 1), c=np.clip(depth[ii], 0, 100), cmap=cm, vmin=0, vmax=100, edgecolors='white', linewidths=0.3, alpha=0.9)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
fig.savefig(f'map_drowned_{SUF}.png', dpi=100); plt.close(fig)
