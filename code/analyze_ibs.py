"""Compare hand-drawn ('box') and LOVECLIM ('ibs') climates on the paired Monte Carlo (mc_ibs.pkl),
and benchmark site arrival times against Timmermann & Friedrich (2016) scenario A/B/C.
Run from the work folder:  python ../code/analyze_ibs.py      -> summary_ibs.json, fig_ibs.png"""
import json, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import model as M, climate_ibs as CI, mc as MC

df = pd.read_pickle('mc_ibs.pkl')
S, CL = ['coastal', 'terrestrial', 'mixed'], ['box', 'ibs']
q = lambda x: [round(float(v), 2) for v in np.nanpercentile(x, [10, 50, 90])] if len(x) else None
out = {'n_runs': {sc: int((df.scen == sc).sum()) for sc in df.scen.unique()}}

# ---------------------------------------------------------------- 1. outcomes by climate
for sc in sorted(df.scen.unique()):
    D = df[df.scen == sc]; o = {}
    for cl in CL:
        r = {}
        for s in S:
            ts = D[f'{cl}|{s}|Sahul']
            ok = np.isfinite(ts)
            r[s] = dict(reach_sahul=round(float(ok.mean()), 3),
                        transit_ka=q(((D.t0 - ts) / 1000)[ok]),
                        route_sunk=q(D.loc[ok, f'{cl}|{s}|route_sunk']))
        c = D[f'{cl}|coastal|Sahul'].replace(-np.inf, 0); t = D[f'{cl}|terrestrial|Sahul'].replace(-np.inf, 0)
        win = (c > t)
        r['coastal_beats_terrestrial'] = round(float(win.mean()), 3)
        qs = pd.qcut(D.c_adv, 4, labels=False)
        r['coastal_win_by_c_adv_quartile'] = [round(float(win[qs == k].mean()), 3) for k in range(4)]
        qb = pd.qcut(D.biome_pen, 4, labels=False)
        r['coastal_win_by_biome_pen_quartile'] = [round(float(win[qb == k].mean()), 3) for k in range(4)]
        o[cl] = r
    cb = (D['box|coastal|Sahul'].replace(-np.inf, 0) > D['box|terrestrial|Sahul'].replace(-np.inf, 0))
    ci = (D['ibs|coastal|Sahul'].replace(-np.inf, 0) > D['ibs|terrestrial|Sahul'].replace(-np.inf, 0))
    o['same_winner_both_climates'] = round(float((cb == ci).mean()), 3)
    out[f'outcomes_{sc}'] = o

# ---------------------------------------------------------------- 2. evidence filter (rejection ABC), long chronology
L = df[df.scen == 'L']
if len(L):
    abc = {}
    for cl in CL:
        ev = {s: (L[f'{cl}|{s}|Sulawesi (Maros / Muna)'] >= 67800) & (L[f'{cl}|{s}|Nefud palaeolakes (N Arabia)'] >= 70000)
                 & (L[f'{cl}|{s}|Lida Ajer (Sumatra)'] >= 63000) for s in S}
        tot = sum(int(ev[s].sum()) for s in S)
        acc = pd.concat([L[ev[s]] for s in S]) if tot else L.iloc[:0]
        a = dict(accepted={s: int(ev[s].sum()) for s in S},
                 post_strategy={s: (round(ev[s].sum() / tot, 3) if tot else None) for s in S},
                 params={prm: {'prior': q(L[prm] * sc_), 'posterior': q(acc[prm] * sc_)}
                         for prm, sc_ in [('c_adv', 1), ('v_max', 1), ('t0', 1e-3), ('fw_floor', 1), ('biome_pen', 1)]})
        abc[cl] = a
    out['abc_L'] = abc

# ---------------------------------------------------------------- 3. benchmark vs Timmermann & Friedrich arrival times
try:
    from netCDF4 import Dataset
    ds = Dataset(os.path.join(M.ROOT, 'data', 'external', 'ibs_hds', 'arrival-time.nc'))
    lat = np.asarray(ds['lat'][:], float); lon = np.asarray(ds['lon'][:], float)
    pts = {'S Arabia coast (Yemen)': (13.0, 45.0), 'Jebel Faya (UAE)': (25.1, 55.8), 'Nefud palaeolakes (N Arabia)': (28.0, 40.5),
           'Levant (Skhul/Manot)': (32.6, 35.4), 'Jwalapuram (S India)': (15.3, 78.1), 'Lida Ajer (Sumatra)': (-0.25, 100.3),
           'Sulawesi (Maros / Muna)': (-5.0, 119.7), 'Madjedbebe (N Australia)': (-12.5, 132.9)}
    bench = {}
    for v in ('arrivaltime_a', 'arrivaltime_b', 'arrivaltime_c'):
        if v not in ds.variables: continue
        a = np.ma.filled(np.ma.masked_invalid(np.squeeze(np.asarray(ds[v][:], float))), np.nan)
        a = np.where((a <= 0) | (np.abs(a) > 1e10), np.nan, a)
        LA, LO = np.meshgrid(lat, lon, indexing='ij'); good = np.isfinite(a)
        r = {}
        for name, (la, lo) in pts.items():
            d = (LA - la) ** 2 + ((LO - lo) * np.cos(np.radians(la))) ** 2
            d = np.where(good, d, np.inf); i = np.unravel_index(np.argmin(d), d.shape)
            r[name] = round(float(a[i]), 1) if np.isfinite(d[i]) and d[i] < 9 else None
        bench[v] = r
    ours = {}
    for cl in CL:
        for s in S:
            ours[f'{cl}|{s}'] = {name: (round(float(np.median(x[np.isfinite(x)]) / 1000), 1) if np.isfinite(x).mean() >= 0.5 else None)
                                 for name in pts for x in [L[f'{cl}|{s}|{name}'].values]}
    out['benchmark'] = {'timmermann_arrivaltime_ka': bench,
                        'note_arrivaltime': "Timmermann 'arrivaltime' = time since LAST CONTINUOUS settlement, not first "
                                            "arrival; not comparable to our first-arrival times. Kept for reference only."}

    # proper benchmark: FIRST arrival from the density snapshots (oldest 1-kyr snapshot with density > threshold)
    dd = Dataset(os.path.join(M.ROOT, 'data', 'external', 'ibs_hds', 'human_density.nc'))
    dlat = np.asarray(dd['lat'][:], float); dlon = np.asarray(dd['lon'][:], float); dmt = np.asarray(dd['mtime'][:], float)
    DLA, DLO = np.meshgrid(dlat, dlon, indexing='ij')
    first = {}
    for v in ('human_density_a', 'human_density_b'):
        den = np.ma.filled(np.ma.masked_invalid(np.asarray(dd[v][:], float)), np.nan)     # (time, lat, lon)
        den = np.where(np.abs(den) > 1e10, np.nan, den)
        land_d = np.isfinite(den).any(0)
        for thr in (0.01, 0.1):
            r = {}
            for name, (la, lo) in pts.items():
                dist = (DLA - la) ** 2 + ((DLO - lo) * np.cos(np.radians(la))) ** 2
                dist = np.where(land_d, dist, np.inf); i = np.unravel_index(np.argmin(dist), dist.shape)
                ser = den[:, i[0], i[1]]; occ = np.where(np.nan_to_num(ser) > thr)[0]
                r[name] = float(dmt[occ].max()) if len(occ) and dist[i] < 9 else None
            first[f'{v}|thr{thr}'] = r
    out['benchmark']['timmermann_first_arrival_ka'] = first

    # compare TIMING RELATIVE TO THE FIRST SITE (S Arabia coast) so different exit dates don't matter:
    # lag_ka(site) = arrival(S Arabia) - arrival(site), ours per run (median over runs reaching both), theirs from density
    ref = 'S Arabia coast (Yemen)'; comp = {}
    for cl in CL:
        for s in S:
            a0 = L[f'{cl}|{s}|{ref}'].values
            ours_lag = {}
            for name in pts:
                if name == ref: continue
                a1 = L[f'{cl}|{s}|{name}'].values; ok = np.isfinite(a0) & np.isfinite(a1)
                ours_lag[name] = (round(float(np.median((a0[ok] - a1[ok]) / 1000)), 1), round(float(ok.mean()), 2)) if ok.mean() >= 0.25 else (None, round(float(ok.mean()), 2))
            row = {'ours_lag_ka_and_share': ours_lag}
            for key, r in first.items():
                if r.get(ref) is None: continue
                th = {n: r[ref] - r[n] for n in pts if n != ref and r.get(n) is not None}
                pairs = [(th[n], ours_lag[n][0]) for n in th if ours_lag[n][0] is not None]
                if len(pairs) >= 4:
                    x, y = np.array(pairs).T
                    row[key] = dict(n=len(pairs), spearman=round(float(pd.Series(x).corr(pd.Series(y), method='spearman')), 2),
                                    median_abs_diff_ka=round(float(np.median(np.abs(y - x))), 1),
                                    theirs_lag_ka={n: round(th[n], 1) for n in th})
            comp[f'{cl}|{s}'] = row
    out['benchmark']['lag_vs_first_arrival'] = comp
    out['benchmark']['note_lag'] = ("lag = years after reaching the S Arabian coast; ours = median over L runs that reach both "
                                    "(share = fraction of runs reaching both; None if <25%); theirs from density first-arrival. "
                                    "Not independent validation: both models share LOVECLIM-family climate assumptions.")
except Exception as e:
    import traceback
    out['benchmark'] = {'error': repr(e), 'trace': traceback.format_exc()[-800:]}

# ---------------------------------------------------------------- 4. climate diagnostic: Arabian humidity in LOVECLIM
C = CI.load(); land = ~M.water_mask(-60.0)
belt = land & (M.LAT2 > 16) & (M.LAT2 < 30) & (M.LON2 > 35) & (M.LON2 < 58)
out['loveclim_arabia_desert_fr_pct'] = {int(k): round(float(C['desert_fr'][i][belt].mean()), 1)
                                        for i, k in enumerate(C['ka']) if 35 <= k <= 100}

json.dump(out, open('summary_ibs.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

# ---------------------------------------------------------------- figure
if not len(L):
    raise SystemExit('no long-chronology runs; skipping figure')
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
col = {'box': '#9aa5b1', 'ibs': '#2a78d6'}
o = out.get('outcomes_L', {})
for j, cl in enumerate(CL):
    if cl in o:
        ax[0].bar(np.arange(4) + (j - 0.5) * 0.38, o[cl]['coastal_win_by_c_adv_quartile'], 0.38, color=col[cl], label=cl)
ax[0].set_xticks(range(4), ['Q1 poor shore', 'Q2', 'Q3', 'Q4 rich shore']); ax[0].set_ylabel('coastal beats megafauna'); ax[0].set_ylim(0, 1)
ax[0].legend(frameon=False); ax[0].set_title('Coastal wins by coastal productivity (L)')
data, labels, colors = [], [], []
for s in S:
    for cl in CL:
        ts = L[f'{cl}|{s}|Sahul']; ok = np.isfinite(ts)
        data.append(((L.t0 - ts) / 1000)[ok].values); labels.append(f'{s[:5]}\n{cl}'); colors.append(col[cl])
b = ax[1].boxplot(data, patch_artist=True, showfliers=False, whis=(10, 90))
for patch, c in zip(b['boxes'], colors): patch.set_facecolor(c)
ax[1].set_xticks(range(1, len(labels) + 1), labels, fontsize=8); ax[1].set_ylabel('kyr Horn -> Sahul'); ax[1].set_title('Transit time (L)')
d = out['loveclim_arabia_desert_fr_pct']; ks = sorted(d)
ax[2].plot(ks, [d[k] for k in ks], color='#2a78d6'); ax[2].invert_xaxis()
for a_, b_ in [(105, 95), (86, 70), (60, 52)]: ax[2].axvspan(a_, b_, color='#d8e8d0', alpha=0.6, lw=0)
ax[2].set_xlabel('ka'); ax[2].set_ylabel('Arabia desert fraction %'); ax[2].set_title('Arabia: LOVECLIM vs hand-set wet pulses')
fig.tight_layout(); fig.savefig('fig_ibs.png', dpi=110)
print('wrote summary_ibs.json, fig_ibs.png')
