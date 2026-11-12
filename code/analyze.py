import numpy as np, pandas as pd
df=pd.read_pickle('mc.pkl'); S=['coastal','terrestrial','mixed']
def ok(col,by): return (df[col]>=by)
for scen in 'LS':
    d=df[df.scen==scen]; print('\n=== scenario',scen,len(d))
    for s in S:
        sah=d[f'{s}|Sahul']; sul=d[f'{s}|Sulawesi (Maros / Muna)']
        tr=(d.t0-sah)[np.isfinite(sah)]
        row=dict(reach_sahul_ever=np.isfinite(sah).mean(), transit_med=tr.median(), transit_q25=tr.quantile(.25), transit_q75=tr.quantile(.75),
                 sunk_med=d[f'{s}|route_sunk'].median())
        if scen=='L':
            row.update(sul_by_678=(sul>=67800).mean(), sahul_by_65=(sah>=65000).mean(),
                       lida_by_63=(d[f'{s}|Lida Ajer (Sumatra)']>=63000).mean(),
                       nefud_70_85=(d[f'{s}|Nefud palaeolakes (N Arabia)']>=70000).mean())
        else:
            row.update(sahul_by_45=(sah>=45000).mean(), sahul_by_43=(sah>=43000).mean(),
                       sul_by_45=(sul>=45000).mean())
        print(s,{k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items()})
