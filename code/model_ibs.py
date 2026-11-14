"""
Same dispersal model, but vegetation comes from the LOVECLIM transient climate (Timmermann & Friedrich 2016)
instead of hand-drawn biome boxes + one Saharo-Arabian humidity index + one monsoon index.

Both uses of 'biome' are now climate-driven, per cell and per 1-kyr slice:

  1. big-game habitat quality
      sT = 0.02 + 0.88 * q_npp * q_dry * q_cold * (1 - forest_pen * forest)
        q_npp  = clip(NPP / (npp_ref * NPP90), 0, 1)              more plant growth -> more game
        q_dry  = (1 - desert_fr/100) ** des_exp                    desert fraction kills habitat
        q_cold = clip((T - t_min) / 6, 0, 1)                       cold limit (matters only in the far north)
        forest = clip((NPP - npp_forest*NPP90) / (0.3*NPP90), 0,1) closed forest is poor for big game
  2. re-tooling biome classes (a population pays biome_pen years when it moves into a different class)
        0 hyper-arid   desert_fr >= 60 %
        1 semi-arid    desert_fr >= 25 %  or  NPP < 0.25 * NPP90
        3 closed forest NPP >= npp_forest * NPP90 and T >= 18 °C
        2 savanna / open woodland otherwise
      + 10 inside the Sahul faunal realm (east of Wallacea), so entering Sahul always means new fauna.
      The class map changes through time with the climate (the solver reads the class of the slice in which a
      move happens).

Still hand-placed, as in model.py, because a 1° climate model cannot resolve them: permanent-water refuges along
the Nile, Tigris-Euphrates, Indus and the exposed Gulf basin, and the altitude penalty. Sea level, coastlines,
sea crossings and the solver logic are unchanged.

coastal fresh water (dry-phase springs): fw = fw_floor + (1 - fw_floor) * clip((1 - desert_fr/100) / 0.7, 0, 1)
climate chronology uncertainty: the whole LOVECLIM record is shifted by cl_shift (±2 kyr)

New sampled parameters (sample_extra) are drawn from a separate RNG so the base draw is identical to
the original Monte Carlo for the same run index. run(..., biome_mode='static') reproduces the earlier
LOVECLIM variant that kept the hand-drawn re-tooling classes (for comparison).
"""
import heapq
import numpy as np
import numba as nb
from scipy import ndimage
import model as M
import climate_ibs as CI

SAHUL_REALM = ((M.LON2 >= 131.5) & (M.LAT2 <= 0)) | ((M.LAT2 <= -11.2) & (M.LON2 >= 121))


def sample_extra(rng):
    return dict(
        npp_ref=rng.uniform(0.5, 1.0),       # NPP (as fraction of NPP90) that counts as 'best' big-game habitat
        des_exp=rng.uniform(0.5, 2.0),       # how strongly desert fraction suppresses habitat
        t_min=rng.uniform(-4.0, 4.0),        # mean-annual °C below which habitat fails
        npp_forest=rng.uniform(1.0, 1.4),    # NPP (x NPP90) above which closed forest starts
        forest_pen=rng.uniform(0.3, 0.7),    # max quality loss in closed forest for big-game hunters
        cl_shift=rng.uniform(-2000, 2000),   # LOVECLIM chronology shift (yr)
    )


def biome_classes(npp, des, tmp, p, npp90):
    c = np.full(npp.shape, 2, np.int8)
    c[(des >= 25) | (npp < 0.25 * npp90)] = 1
    c[des >= 60] = 0
    c[(npp >= p['npp_forest'] * npp90) & (tmp >= 18)] = 3
    return (c + 10 * SAHUL_REALM).astype(np.int8)


def speed_surfaces(p, strategy, noise):
    """returns speed[K,N] (km/yr), coast[K,N] and biome class[K,N] (time-varying)."""
    C = CI.load(); npp90 = float(C['NPP90'])
    sp = np.zeros((M.K, M.N), np.float32); coastk = np.zeros((M.K, M.N), np.bool_)
    biok = np.zeros((M.K, M.N), np.int8)
    for k in range(M.K):
        bp = M.T_K[k]
        sl = M.sea_level(bp, p['sl_off'], p['sl_jit'])
        water = M.water_mask(sl)
        land = ~water
        coast = land & ndimage.binary_dilation(water & ~M.LAKE, structure=np.ones((3, 3)))
        npp, des, tmp = CI.at(bp, p['cl_shift'])
        q_npp = np.clip(npp / (p['npp_ref'] * npp90), 0, 1)
        q_dry = np.clip(1 - des / 100.0, 0, 1) ** p['des_exp']
        q_cold = np.clip((tmp - p['t_min']) / 6.0, 0, 1)
        forest = np.clip((npp - p['npp_forest'] * npp90) / (0.3 * npp90), 0, 1)
        sT = 0.02 + 0.88 * q_npp * q_dry * q_cold * (1 - p['forest_pen'] * forest)
        sT = ndimage.gaussian_filter(sT, 1.0)
        sT = np.where(M.RIVER, np.maximum(sT, 0.6), sT)
        sT = np.where(M.GULF_BASIN, np.maximum(sT, 0.45), sT)
        sT = sT * M.ALT_PEN
        fw = p['fw_floor'] + (1 - p['fw_floor']) * np.clip((1 - des / 100.0) / 0.7, 0, 1)
        fw = ndimage.gaussian_filter(fw, 1.0)
        sC = p['c_adv'] * fw
        if strategy == 'coastal':
            s = np.where(coast, sC, p['inland_frac'] * sT)
        elif strategy == 'terrestrial':
            s = sT
        else:
            # c_mix = 0 (default) reproduces the published LOVECLIM runs; c_mix > 0 = cost of flexibility (review item 1)
            s = (1.0 - p.get('c_mix', 0.0)) * np.where(coast, np.maximum(sC, sT), sT)
        s = s * noise
        s = np.where(land & (s >= 0.12), s, 0.0)
        sp[k] = (p['v_max'] * s).ravel(); coastk[k] = coast.ravel()
        biok[k] = biome_classes(npp, des, tmp, p, npp90).ravel()
    return sp, coastk, biok


# time-dependent Dijkstra, identical to model._dijkstra except biome is [K, N] and read at the move's slice
@nb.njit(cache=True)
def _dijkstra_tv(speed, coastk, biome, scode, pen, nbr_d, offs, ny, nx, seeds, tau0, slbin, cptr, cidx, cw, lam, w0, dmax, dt, tau_max):
    Kk, Nn = speed.shape
    arr = np.full(Nn, np.inf); pred = np.full(Nn, -1, np.int64); via_sea = np.zeros(Nn, np.bool_)
    heap = [(0.0, np.int64(0))]; heap.pop()
    for s in seeds:
        k0 = int(tau0 // dt)
        if speed[k0, s] > 0:
            arr[s] = tau0; heap.append((tau0, s))
    heapq.heapify(heap)
    while len(heap) > 0:
        t, u = heapq.heappop(heap)
        if t > arr[u] or t > tau_max: continue
        k0 = int(t // dt)
        if k0 >= Kk: continue
        r = u // nx; c = u % nx
        for e in range(8):
            r2 = r + offs[e, 0]; c2 = c + offs[e, 1]
            if r2 < 0 or r2 >= ny or c2 < 0 or c2 >= nx: continue
            v = r2 * nx + c2
            k = k0
            while k < Kk:
                su = speed[k, u]
                if su <= 0: break
                sv = speed[k, v]
                if sv > 0:
                    dep = max(t, k * dt)
                    na = dep + nbr_d[r, e] / (0.5 * (su + sv))
                    if biome[k, u] != biome[k, v] and not (scode != 1 and coastk[k, u] and coastk[k, v]):
                        na += pen
                    if na < arr[v]:
                        arr[v] = na; pred[v] = u; via_sea[v] = False
                        heapq.heappush(heap, (na, v))
                    break
                k += 1
        b = slbin[k0]
        if speed[k0, u] > 0:
            for q in range(cptr[b, u], cptr[b, u + 1]):
                v = cidx[q]; w = cw[q]
                if w > dmax or speed[k0, v] <= 0: continue
                na = t + w0 * np.exp(w / lam)
                if na < arr[v]:
                    arr[v] = na; pred[v] = u; via_sea[v] = True
                    heapq.heappush(heap, (na, v))
    return arr, pred, via_sea


def run(p, strategy, seeds, t0_bp, noise, X, tau_max_bp=M.T_END, biome_mode='loveclim'):
    sp, coastk, biok = speed_surfaces(p, strategy, noise)
    if biome_mode == 'static':
        biok = np.broadcast_to(M.BIOME.ravel().astype(np.int8), biok.shape).copy()
    scode = {'coastal': 0, 'terrestrial': 1, 'mixed': 2}[strategy]
    slb = np.array([np.argmin(np.abs(M.SL_BINS - M.sea_level(M.T_K[k], p['sl_off'], p['sl_jit']))) for k in range(M.K)])
    arr, pred, via = _dijkstra_tv(sp, coastk, biok, scode, float(p['biome_pen']), M.NBR_D, M.OFFS, M.NY, M.NX,
                                  np.asarray(seeds, np.int64), float(M.bp_to_tau(t0_bp)), slb, X[0], X[1], X[2],
                                  p['lam'], p['w0'], p['dmax'], float(M.DT), float(M.bp_to_tau(tau_max_bp)))
    return M.tau_to_bp(arr), pred, via


def run_static_biome(p, strategy, seeds, t0_bp, noise, X, tau_max_bp=M.T_END):
    return run(p, strategy, seeds, t0_bp, noise, X, tau_max_bp, biome_mode='static')


def habitat_snapshot(p, bp):
    """Big-game habitat quality (before noise) at one time, for maps/diagnostics."""
    C = CI.load(); npp90 = float(C['NPP90'])
    npp, des, tmp = CI.at(bp, p.get('cl_shift', 0.0))
    q = (np.clip(npp / (p['npp_ref'] * npp90), 0, 1) * np.clip(1 - des / 100, 0, 1) ** p['des_exp']
         * np.clip((tmp - p['t_min']) / 6, 0, 1)
         * (1 - p['forest_pen'] * np.clip((npp - p['npp_forest'] * npp90) / (0.3 * npp90), 0, 1)))
    return 0.02 + 0.88 * q


def biome_snapshot(p, bp):
    C = CI.load(); npp, des, tmp = CI.at(bp, p.get('cl_shift', 0.0))
    return biome_classes(npp, des, tmp, p, float(C['NPP90']))
