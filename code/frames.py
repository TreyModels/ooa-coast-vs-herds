"""Re-run a subset of long-chronology draws keeping full arrival-time grids, for occupancy-through-time."""
import numpy as np, model as M, mc, os
from multiprocessing import Pool
def one(i):
    rng=np.random.default_rng(i); p=M.sample_params(rng); noise=M.noise_field(rng); t0=rng.uniform(72_000,85_000)
    return [M.run(p,s,mc.HORN,t0,noise,mc.X)[0].astype(np.float32) for s in mc.STRATS]
if __name__=='__main__':
    n=150
    with Pool(os.cpu_count()) as pool: res=pool.map(one,range(n),chunksize=5)
    A=np.array(res)  # n,3,N
    np.save('arr_L.npy',A); print(A.shape)
