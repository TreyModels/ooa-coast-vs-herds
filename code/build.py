import numpy as np, time, model as M
t=time.time(); X=M.build_crossings(); print('edges',len(X[1]), time.time()-t)
np.savez('cross.npz', ptr=X[0], idx=X[1], w=X[2])
# sanity: key straits at -80 m
w=M.water_mask(-80.0)
def gap(name, la0,lo0,la1,lo1):
    ns=200; la=np.linspace(la0,la1,ns); lo=np.linspace(lo0,lo1,ns)
    ii=[int(np.argmin(abs(M.LAT-a))) for a in la]; jj=[int(np.argmin(abs(M.LON-b))) for b in lo]
    seg=w[ii,jj]; km=np.hypot((la1-la0)*111,(lo1-lo0)*111*np.cos(np.radians(la0)))
    print(f"{name}: water fraction {seg.mean():.2f} of {km:.0f} km line")
gap('Bab el Mandeb', 12.9,43.0, 12.7,43.6)
gap('Sinai land bridge', 30.5,32.3, 30.5,34.0)
gap('Gulf floor (Qatar->Iran)', 25.5,52, 27.5,52)
gap('Hormuz', 26.0,56.3, 26.9,56.3)
gap('Sumatra-Java (Sunda shelf)', -5.5,105.6,-6.2,105.9)
gap('Bali-Lombok', -8.5,115.6, -8.5,116.1)
gap('Borneo-Sulawesi (Makassar)', -1.0,117.3, -1.0,119.5)
gap('Timor-Sahul', -10.3,125.5, -11.8,127.5)
print('land cells at -80 vs 0:', (~w).sum(), (~M.water_mask(0.0)).sum())
