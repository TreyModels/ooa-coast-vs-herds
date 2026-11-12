import numpy as np, pandas as pd, json, model as M, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
df=pd.read_pickle('mc.pkl'); D=np.load('density.npz'); n=int(D['n'])
S=['coastal','terrestrial','mixed']; COL={'coastal':'#2a78d6','terrestrial':'#eb6834','mixed':'#1baf7a'}
r=(M.LAT<=40)&(M.LAT>=-15); c=(M.LON>=30)&(M.LON<=150)
ext=[M.LON[c][0],M.LON[c][-1],M.LAT[r][-1],M.LAT[r][0]]
land_now=~M.water_mask(0.0); land_lo=~M.water_mask(-80.0)
sites={'Jebel Faya':(25.1,55.8),'Nefud lakes':(28.0,40.5),'Skhul/Manot':(32.6,35.4),'Jwalapuram':(15.3,78.1),
 'Lida Ajer':(-0.25,100.3),'Sulawesi art':(-4.9,121.5),'Madjedbebe':(-12.5,132.9),'Bab-el-Mandeb':(12.6,43.4)}
def base(ax):
    rgb=np.ones(land_now.shape+(3,))*np.array([0.93,0.95,0.97])
    rgb[land_lo]=[0.80,0.78,0.73]; rgb[land_lo&~land_now]=[0.66,0.62,0.55]; rgb[land_now]=[0.89,0.87,0.83]
    ax.imshow(rgb[np.ix_(r,c)],extent=ext,interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
def marks(ax):
    for k,(la,lo) in sites.items():
        ax.plot(lo,la,'o',ms=4,mfc='#0b0b0b',mec='white',mew=0.8)
        ax.annotate(k,(lo,la),xytext=(4,3),textcoords='offset points',fontsize=7,color='#0b0b0b')
# map 1: exposed shelf
fig,ax=plt.subplots(figsize=(9,4.6)); base(ax); marks(ax); plt.tight_layout(pad=0.2); plt.savefig('map_shelf.png',dpi=110); plt.close()
# map 2/3: route density
for s in S:
    dens=D[f'L|{s}'].reshape(M.NY,M.NX)/n
    fig,ax=plt.subplots(figsize=(9,4.6)); base(ax)
    a=np.ma.masked_less(dens,0.02)[np.ix_(r,c)]
    from matplotlib.colors import LinearSegmentedColormap
    cm=LinearSegmentedColormap.from_list('x',['#ffffff00',COL[s]])
    cm=ListedColormap(np.c_[np.tile(matplotlib.colors.to_rgb(COL[s]),(256,1)),np.linspace(0.25,1,256)])
    ax.imshow(a,extent=ext,cmap=cm,vmin=0.02,vmax=0.5,interpolation='nearest'); marks(ax)
    plt.tight_layout(pad=0.2); plt.savefig(f'map_route_{s}.png',dpi=110); plt.close()
# summaries
L=df[df.scen=='L']; Sd=df[df.scen=='S']
out={'n':n}
T=np.arange(100,34,-1)*1000.0
sl=np.array([M.sea_level(T,o,j) for o,j in zip(L.sl_off,L.sl_jit)])
H=np.array([M.climate(T,row)[0] for _,row in L.iterrows()])
out['timeline']={'ka':(T/1000).tolist(),'sl':[np.percentile(sl,q,axis=0).round(1).tolist() for q in (10,50,90)],
                 'H':[np.percentile(H,q,axis=0).round(3).tolist() for q in (10,50,90)]}
def q(x,ps=(10,25,50,75,90)): x=x[np.isfinite(x)]; return [round(float(np.percentile(x,p))/1000,2) for p in ps]
out['transit']={s:q((L.t0-L[f'{s}|Sahul']).values) for s in S}
out['sunk']={s:q(L[f'{s}|route_sunk'].values*1000) for s in S}
ev={}
for s in S:
    sul=L[f'{s}|Sulawesi (Maros / Muna)']>=67800; nef=L[f'{s}|Nefud palaeolakes (N Arabia)']>=70000
    ev[s]=dict(sul=sul.mean(),sah65=(L[f'{s}|Sahul']>=65000).mean(),lida=(L[f'{s}|Lida Ajer (Sumatra)']>=63000).mean(),
       nef=nef.mean(),joint=(sul&nef).mean(),S45=(Sd[f'{s}|Sahul']>=45000).mean(),S43=(Sd[f'{s}|Sahul']>=43000).mean(),
       Sever=np.isfinite(Sd[f'{s}|Sahul']).mean())
out['evidence']={s:{k:round(float(v),3) for k,v in e.items()} for s,e in ev.items()}
c_=L['coastal|Sahul'].replace(-np.inf,0); t_=L['terrestrial|Sahul'].replace(-np.inf,0); win=(c_>t_)
out['win_overall']=round(float(win.mean()),3); out['delta_med']=round(float((c_-t_).median()))
def pdp(prm,labels):
    qq=pd.qcut(L[prm],4); g=win.groupby(qq,observed=True).mean()
    return [{'label':lab,'range':f"{iv.left:.2f}–{iv.right:.2f}",'win':round(float(v),3)} for lab,(iv,v) in zip(labels,g.items())]
out['pdp']={'c_adv':pdp('c_adv',['Q1','Q2','Q3','Q4']),'biome_pen':pdp('biome_pen',['Q1','Q2','Q3','Q4']),'fw_floor':pdp('fw_floor',['Q1','Q2','Q3','Q4'])}
fav=(L.c_adv>1.0)&(L.fw_floor>0.35)&(L.biome_pen>700); unf=(L.c_adv<0.9)&(L.biome_pen<500)
out['corners']={'fav':[int(fav.sum()),round(float(win[fav].mean()),2)],'unf':[int(unf.sum()),round(float(win[unf].mean()),2)]}
json.dump(out,open('summary.json','w'))
print(json.dumps({k:out[k] for k in ['transit','sunk','evidence','win_overall','delta_med','corners']},indent=0))
