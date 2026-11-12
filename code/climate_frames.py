import numpy as np, model as M, matplotlib, io, base64, json
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from scipy import ndimage
from PIL import Image
A=np.load('arr_L.npy')  # bp arrival, -inf/inf if never
S=['coastal','terrestrial','mixed']; COL={'coastal':'#2a78d6','terrestrial':'#eb6834','mixed':'#1baf7a'}
r=(M.LAT<=40)&(M.LAT>=-15); c=(M.LON>=30)&(M.LON<=150)
ext=[M.LON[c][0],M.LON[c][-1],M.LAT[r][-1],M.LAT[r][0]]
# representative (mid-range) climate parameters
p=dict(sl_off=0,sl_jit=0,H_bg=0.07,h5c_a=103.5,h5c_b=96.5,h5c_amp=0.75,h5a_a=84,h5a_b=71.5,h5a_amp=0.75,
       h3_a=58.5,h3_b=53.5,h3_amp=0.3,H_mis4=0.06,M_mis4=0.4,M_mis3=0.68)
hab=LinearSegmentedColormap.from_list('hab',['#f6efe1','#7a5a2e'])
land_now=~M.water_mask(0.0)
times=np.arange(85000,44000,-2500)
frames={s:[] for s in S}; meta=[]
for t in times:
    sl=float(M.sea_level(t)); H,Mi=M.climate(np.array([float(t)]),p); h,m=H[0],Mi[0]
    water=M.water_mask(sl); land=~water; glacial=np.clip(-sl/100,0,1)
    B=M.BIOME
    sT=np.select([B==0,B==1,B==2,B==3,B==4],[0.02+0.85*h,0.3+0.45*h,0.45+0.4*m,0.25+0.35*glacial,0.8])
    sT=ndimage.gaussian_filter(sT,1.5); sT=np.where(M.RIVER,np.maximum(sT,0.6),sT); sT=np.where(M.GULF_BASIN,np.maximum(sT,0.45),sT)*M.ALT_PEN
    rgba=hab(np.clip(sT/0.9,0,1)); rgba[water]=to_rgb('#dbe7f1')+(1,)
    shelf=land&~land_now
    meta.append(dict(ka=t/1000,sl=round(sl),H=round(float(h),2)))
    for k,s in enumerate(S):
        occ=(A[:,k,:]>=t).mean(0).reshape(M.NY,M.NX)
        fig=plt.figure(figsize=(9,4.12)); ax=fig.add_axes([0,0,1,1]); ax.set_axis_off()
        ax.imshow(rgba[np.ix_(r,c)],extent=ext,interpolation='nearest')
        # present-day coastline on the exposed shelf, for reference
        ax.contour(M.LON[c],M.LAT[r],land_now[np.ix_(r,c)].astype(float),[0.5],colors=['#8a8070'],linewidths=0.5,linestyles=':')
        o=occ[np.ix_(r,c)]
        col=np.array(to_rgb(COL[s]))
        import matplotlib.patheffects as pe
        o=ndimage.gaussian_filter(o,1.2)
        if (o>=0.5).any():
            cs=ax.contour(M.LON[c],M.LAT[r],o,[0.5],colors=[COL[s]],linewidths=2.2)
            for coll in ([cs] if not hasattr(cs,'collections') else cs.collections): coll.set_path_effects([pe.Stroke(linewidth=4.2,foreground='white'),pe.Normal()])
        buf=io.BytesIO(); fig.savefig(buf,dpi=100,format='png'); plt.close(fig)
        im=Image.open(buf).convert('RGB'); b2=io.BytesIO(); im.save(b2,'WEBP',quality=78)
        frames[s].append('data:image/webp;base64,'+base64.b64encode(b2.getvalue()).decode())
json.dump({'meta':meta,'frames':frames},open('frames.json','w'))
print(sum(len(f) for v in frames.values() for f in v)/1e6,'MB')
Image.open(io.BytesIO(base64.b64decode(frames['coastal'][8].split(',')[1]))).save('frame_check.png')

