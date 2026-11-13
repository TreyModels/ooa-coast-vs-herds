import numpy as np, pandas as pd, json, erectus as E
df = pd.read_pickle('erectus.pkl'); out = {'n': len(df)}
tab = {}
for s in E.STRATS:
    for c in E.CROSS:
        k = f'{s}|{c}'; r = {}
        for name, (la, lo, age, rad) in E.PRESENT.items():
            r[name] = round(float((df[f'{k}|{name}'] >= age).mean()), 3)
        for name in E.ABSENT:
            r['NOT ' + name] = round(float((df[f'{k}|NOT {name}'] == -np.inf).mean()), 3)
        pres = np.all([df[f'{k}|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items()], axis=0)
        mainland = np.all([df[f'{k}|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items() if 'Sulawesi' not in n and 'Flores' not in n and 'Luzon' not in n], axis=0)
        islands = np.all([df[f'{k}|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items() if any(x in n for x in ['Sulawesi','Flores','Luzon'])], axis=0)
        absent = np.all([df[f'{k}|NOT {n}'] == -np.inf for n in E.ABSENT], axis=0)
        r['_mainland_all'] = round(float(mainland.mean()), 3); r['_islands_all'] = round(float(islands.mean()), 3)
        r['_absent_all'] = round(float(absent.mean()), 3); r['_everything'] = round(float((pres & absent).mean()), 3)
        tab[k] = r
out['table'] = tab
# timing: how fast to the earliest sites relative to start
out['lag_shangchen'] = {s: [round(float(np.percentile((df.t0 - df[f'{s}|walk|Shangchen, Loess Plateau (2.1 Ma)']).replace(np.inf, np.nan).dropna(), q)) / 1000) for q in (10, 50, 90)] for s in E.STRATS}
out['lag_dmanisi'] = {s: [round(float(np.percentile((df.t0 - df[f'{s}|walk|Dmanisi, Georgia (1.85 Ma)']).replace(np.inf, np.nan).dropna(), q)) / 1000) for q in (10, 50, 90)] for s in E.STRATS}
# canoe: what else, and when
canoe_first = {}
for name in E.ABSENT:
    v = df[f'mixed|canoe|NOT {name}']; v = v[np.isfinite(v)]
    canoe_first[name] = {'p': round(float(len(v) / len(df)), 3), 'median_ka': round(float(v.median()) / 1000) if len(v) else None}
out['canoe_first'] = canoe_first
# sea-level model sensitivity for the raft runs
out['sl_sens'] = {m: {k: round(float((g[f'mixed|raft|{k}'] >= a).mean()), 3) for k, (la, lo, a, rr) in E.PRESENT.items() if any(x in k for x in ['Sulawesi','Flores','Luzon'])} for m, g in df.groupby('sl_model')}
out['sl_sens_sahul'] = {m: round(float((g['mixed|raft|NOT Sahul (New Guinea + Australia)'] == -np.inf).mean()), 3) for m, g in df.groupby('sl_model')}
# raft frequency that fits islands but not Sahul
d = df.copy(); ok_is = np.all([d[f'mixed|raft|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items() if any(x in n for x in ['Sulawesi','Flores','Luzon'])], axis=0)
ok_ab = d['mixed|raft|NOT Sahul (New Guinea + Australia)'] == -np.inf
bins = [20e3, 50e3, 100e3, 200e3, 300e3]
out['raft_w0'] = [{'bin': f'{int(a/1000)}–{int(b/1000)}k', 'islands': round(float(ok_is[(d.raft_w0 >= a) & (d.raft_w0 < b)].mean()), 2),
                   'no_sahul': round(float(ok_ab[(d.raft_w0 >= a) & (d.raft_w0 < b)].mean()), 2),
                   'both': round(float((ok_is & ok_ab)[(d.raft_w0 >= a) & (d.raft_w0 < b)].mean()), 2)} for a, b in zip(bins[:-1], bins[1:])]
json.dump(out, open('er.json', 'w'), indent=1)
for k, r in tab.items(): print(k, {kk: v for kk, v in r.items()})
print(json.dumps({k: out[k] for k in ['lag_shangchen', 'lag_dmanisi', 'canoe_first', 'sl_sens', 'sl_sens_sahul', 'raft_w0']}, indent=1))

# combined fit, leaving out Japan (absence there is weak evidence; a glacial land bridge is debated)
ABS = [n for n in E.ABSENT if 'Japan' not in n]
fit = {}
for s in E.STRATS:
    for c in E.CROSS:
        k = f'{s}|{c}'
        pres = np.all([df[f'{k}|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items()], axis=0)
        absent = np.all([df[f'{k}|NOT {n}'] == -np.inf for n in ABS], axis=0)
        fit[k] = dict(present=round(float(pres.mean()), 3), absent=round(float(absent.mean()), 3), both=round(float((pres & absent).mean()), 3))
out['fit'] = fit
json.dump(out, open('er.json', 'w'), indent=1)
print(json.dumps(fit, indent=0))
# best rafting rate given mixed strategy
d = df; k = 'mixed|raft'
pres = np.all([d[f'{k}|{n}'] >= a for n, (la, lo, a, rr) in E.PRESENT.items()], axis=0)
absent = np.all([d[f'{k}|NOT {n}'] == -np.inf for n in ABS], axis=0)
ok = pres & absent
print('raft w0 of fitting runs (k yr):', np.percentile(d.raft_w0[ok] / 1000, [10, 50, 90]).round(0), 'n', ok.sum())
