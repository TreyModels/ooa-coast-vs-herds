"""Review tests on the LOVECLIM climate (run from work/ after mc_ibs.py runs).

Inputs: mc_ibs.pkl (c_mix 0, both climates), mc_ibs_c0.1.pkl, mc_ibs_c0.2.pkl (LOVECLIM only).
Evidence filter as in analyze_ibs.py: Sulawesi by 67.8 ka, Sumatra (Lida Ajer) by 63 ka, Nefud occupied by 70 ka.
'no Nefud' drops the Nefud condition. Writes ../results/loveclim/review_ibs.json.
"""
import json
import numpy as np, pandas as pd

S = ['coastal', 'terrestrial', 'mixed']


def summarise(df, cl):
    out = {}
    for sc in ('L', 'S'):
        D = df[df.scen == sc]
        if not len(D): continue
        sah = pd.DataFrame({s: D[f'{cl}|{s}|Sahul'].replace(-np.inf, np.nan) for s in S})
        r = {'n': len(D),
             'P(mixed earliest at Sahul)': round(float((sah['mixed'] >= sah[['coastal', 'terrestrial']].max(axis=1)).mean()), 3),
             'P(coastal earlier than terrestrial)': round(float((sah['coastal'].fillna(0) > sah['terrestrial'].fillna(0)).mean()), 3),
             'reach Sahul': {s: round(float(sah[s].notna().mean()), 3) for s in S}}
        if sc == 'L':
            for name, use_nef in (('posterior', True), ('posterior_no_nefud', False)):
                ev = {}
                for s in S:
                    e = (D[f'{cl}|{s}|Sulawesi (Maros / Muna)'] >= 67800) & (D[f'{cl}|{s}|Lida Ajer (Sumatra)'] >= 63000)
                    if use_nef: e &= D[f'{cl}|{s}|Nefud palaeolakes (N Arabia)'] >= 70000
                    ev[s] = int(e.sum())
                tot = sum(ev.values())
                r[name] = {s: round(ev[s] / tot, 3) if tot else None for s in S}
                r[name + '_accepted'] = ev
        out[sc] = r
    return out


res = {}
base = pd.read_pickle('mc_ibs.pkl')
res['box c_mix 0'] = summarise(base, 'box')
res['loveclim c_mix 0'] = summarise(base, 'ibs')
for c in ('0.1', '0.2'):
    try: res[f'loveclim c_mix {c}'] = summarise(pd.read_pickle(f'mc_ibs_c{c}.pkl'), 'ibs')
    except FileNotFoundError: pass
json.dump(res, open('../results/loveclim/review_ibs.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
