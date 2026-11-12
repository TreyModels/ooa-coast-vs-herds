import numpy as np, pandas as pd
df=pd.read_pickle('mc.pkl')
d=df[df.scen=='L'].copy()
c=d['coastal|Sahul'].replace(-np.inf,0); t=d['terrestrial|Sahul'].replace(-np.inf,0)
d['delta']=(c-t)  # +ve: coastal arrives earlier (years)
d['cwin']=d.delta>0
print('coastal earlier in', d.cwin.mean().round(3), 'median delta yrs', d.delta.median().round(0))
for prm in ['c_adv','fw_floor','biome_pen','v_max','h5a_amp','inland_frac','w0','lam','t0','H_mis4']:
    q=pd.qcut(d[prm],4)
    g=d.groupby(q,observed=True).agg(cwin=('cwin','mean'),delta=('delta','median'))
    print(prm, ' | '.join(f"{iv.left:.2f}-{iv.right:.2f}: win {r.cwin:.2f} Δ{r.delta:+.0f}" for iv,r in g.iterrows()))
# joint: high coastal productivity & good coastal freshwater & high biome penalty
m=(d.c_adv>1.0)&(d.fw_floor>0.35)&(d.biome_pen>700)
print('favourable corner n',m.sum(),'coastal wins',d[m].cwin.mean().round(2),'median Δ',d[m].delta.median().round(0))
m2=(d.c_adv<0.9)&(d.biome_pen<500)
print('unfavourable corner n',m2.sum(),'coastal wins',d[m2].cwin.mean().round(2),'median Δ',d[m2].delta.median().round(0))
# exit timing: by t0 bins, where t0 in MIS5a humid vs late
for lo,hi in [(72,76),(76,80),(80,85)]:
    s=d[(d.t0>=lo*1000)&(d.t0<hi*1000)]
    print(f't0 {lo}-{hi}ka n={len(s)} coastal win {s.cwin.mean():.2f}  P(Sul by 67.8) coast {(s["coastal|Sulawesi (Maros / Muna)"]>=67800).mean():.2f} terr {(s["terrestrial|Sulawesi (Maros / Muna)"]>=67800).mean():.2f} mixed {(s["mixed|Sulawesi (Maros / Muna)"]>=67800).mean():.2f}')
print('route crossings med', {s:d[f'{s}|route_crossings'].median() for s in ['coastal','terrestrial','mixed']})
dS=df[df.scen=='S']
for s in ['coastal','terrestrial','mixed']:
    x=dS[np.isfinite(dS[f'{s}|Sahul'])]
    print('S',s,'needed v_max for Sahul by 45 ka:', x[x[f'{s}|Sahul']>=45000].v_max.round(1).tolist())
