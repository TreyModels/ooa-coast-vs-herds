"""Gate and summary for the review reruns (run from work/)."""
import json, sys, glob
import numpy as np, pandas as pd
mode = sys.argv[1]
if mode == 'gate':
    new = pd.read_pickle('mc_c0_w0.pkl'); old = pd.read_csv('../results/monte_carlo_runs.csv')
    assert len(new) == len(old), (len(new), len(old))
    worst = 0.0
    for c in old.columns:
        if c not in new.columns: sys.exit(f'GATE FAIL: column {c} missing')
        if old[c].dtype.kind in 'fi':
            a, b = old[c].to_numpy(float), new[c].to_numpy(float)
            same_nan = np.array_equal(np.isnan(a), np.isnan(b))
            d = np.nanmax(np.abs(a - b) / np.maximum(1, np.abs(a))) if np.isfinite(a).any() else 0
            if not same_nan or d > 1e-6: sys.exit(f'GATE FAIL: {c} differs (rel {d:.2e}, nan-match {same_nan})')
            worst = max(worst, d)
        elif not (old[c].astype(str).values == new[c].astype(str).values).all():
            sys.exit(f'GATE FAIL: {c} differs')
    e_old = json.load(open('../results/exp.json')); e_new = json.load(open('exp_c0_w0.json'))
    for k in ('post_strategy', 'abc'):
        if json.dumps(e_old[k], sort_keys=True) != json.dumps(e_new[k], sort_keys=True):
            sys.exit(f'GATE FAIL: exp.json key {k} differs:\n old {e_old[k]}\n new {e_new[k]}')
    print(f'GATE OK: {len(old)} runs reproduce (worst rel diff {worst:.1e}); exp post_strategy/abc identical')
else:
    rows = {}
    for f in sorted(glob.glob('mc_*.pkl')):
        tag = f[3:-4]; df = pd.read_pickle(f)
        sah = {s: df[f'{s}|Sahul'].replace(-np.inf, np.nan) for s in ('coastal', 'terrestrial', 'mixed') if f'{s}|Sahul' in df}
        S = pd.DataFrame(sah)
        r = {'n': len(df),
             'P(mixed earliest at Sahul)': round(float((S['mixed'] >= S[['coastal', 'terrestrial']].max(axis=1)).mean()), 3),
             'P(coastal earlier than terrestrial)': round(float((S['coastal'] > S['terrestrial']).mean()), 3),
             'Sahul arrival ka, median by strategy': {k: round(float(np.nanmedian(v)) / 1000, 1) for k, v in sah.items()}}
        for suf in ('', '_nonefud'):
            try: e = json.load(open(f'exp_{tag}{suf}.json'))
            except FileNotFoundError: continue
            r['abc' + suf] = e.get('abc'); r['post_strategy' + suf] = e.get('post_strategy')
            if 'hotspots' in e and not suf: r['hotspots'] = e['hotspots']
        rows[tag] = r
    json.dump(rows, open('../results/review_summary.json', 'w'), indent=1, default=str)
    print(json.dumps({t: {k: v for k, v in r.items() if k != 'hotspots'} for t, r in rows.items()}, indent=1, default=str))
