"""
Download ICCP/IBS CESM hominin habitat-suitability layers (Timmermann et al. 2022) and CESM 3 Ma climate
(Yun et al. 2023) over OPeNDAP (GrADS-DODS) into local, compressed netCDF files.

The website's 'Data Access' buttons are dynamic and often do not resolve; the DODS server behind them works:
  https://climatedata.ibs.re.kr:9443/dods/public-data/...   (fallback http://...:9090/...)

Chunked by time, resumable (a chunk already written is skipped), retries with back-off, pauses between
requests to stay gentle on a small server. Run from anywhere:
    python code/fetch_iccp.py              # everything below
    python code/fetch_iccp.py habitat      # only the habitat layers
    python code/fetch_iccp.py climate      # only CESM ts / prect
Writes data/external/iccp_cesm/*.nc, SHA256SUMS and fetch_log.txt.

Citations: Timmermann et al. 2022, Nature 604:495-501, doi:10.1038/s41586-022-04600-9 (habitat layers);
Yun et al. 2023, Clim. Past 19:1951-1974, doi:10.5194/cp-19-1951-2023, data doi:10.22741/iccp.20230001
(3 Ma CESM1.2 run). ICCP licence statement: 'CC BY 4.0 ... CC BY-NC' (contradictory; treat as attribution +
noncommercial). Nature reporting summary: 'no restriction on data usage and access'.
"""
import os, sys, time, hashlib, datetime
import numpy as np
from netCDF4 import Dataset

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'data', 'external', 'iccp_cesm')
HOSTS = ['https://climatedata.ibs.re.kr:9443/dods/public-data', 'http://climatedata.ibs.re.kr:9090/dods/public-data']
CHUNK = 100          # time steps per request
PAUSE = 1.0          # seconds between requests

JOBS = {
    'habitat': dict(path='hominin-habitats/hominin-habitats',
                    vars=['erectus', 'ergahabi', 'heidelberg', 'neander', 'sapiens'],
                    lon=(0, 180), lat=(-40, 70),
                    # 2001 steps; index 0 = 2000 ka ... index 2000 = present (1-kyr steps, model time axis)
                    ka=lambda n: 2000.0 - np.arange(n),
                    cite='Timmermann et al. 2022, Nature 604:495-501, doi:10.1038/s41586-022-04600-9',
                    note='Mahalanobis-distance niche suitability; 0 = missing / species absent. '
                         'ICCP page notes the runs were made without CO2 fertilisation.'),
    'climate': dict(path='pcesm-3ma/pcesm-3ma',
                    vars=['ts', 'prect'],
                    lon=None, lat=None,                       # global, 3.75 deg
                    # 3001 steps; index 0 = 3000 ka ... index 3000 = present
                    ka=lambda n: 3000.0 - np.arange(n),
                    cite='Yun et al. 2023, Clim. Past 19:1951-1974, doi:10.5194/cp-19-1951-2023; '
                         'data doi:10.22741/iccp.20230001',
                    note='CESM1.2 3 Ma transient run, annual means, ~3.75 deg. ts in K (surface temp), '
                         'prect in m/s (total precip rate).'),
}


def log(msg):
    line = f'{datetime.datetime.now():%Y-%m-%d %H:%M:%S}  {msg}'
    print(line, flush=True)
    with open(os.path.join(OUT, 'fetch_log.txt'), 'a') as f:
        f.write(line + '\n')


def open_remote(path):
    last = None
    for attempt in range(6):
        for h in HOSTS:
            try:
                return Dataset(f'{h}/{path}'), f'{h}/{path}'
            except Exception as e:
                last = e
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'cannot open {path}: {last!r}')


def fetch_var(job, name, var):
    out = os.path.join(OUT, f'{name}_{var}.nc')
    src, url = open_remote(job['path'])
    lat = np.asarray(src['lat'][:], float); lon = np.asarray(src['lon'][:], float)
    li = np.arange(len(lat)) if job['lat'] is None else np.flatnonzero((lat >= job['lat'][0]) & (lat <= job['lat'][1]))
    lj = np.arange(len(lon)) if job['lon'] is None else np.flatnonzero((lon >= job['lon'][0]) & (lon <= job['lon'][1]))
    nt = src[var].shape[0]
    i0, i1, j0, j1 = li[0], li[-1] + 1, lj[0], lj[-1] + 1
    if not os.path.exists(out):
        d = Dataset(out, 'w')
        d.createDimension('time', nt); d.createDimension('lat', i1 - i0); d.createDimension('lon', j1 - j0)
        v = d.createVariable('ka_bp', 'f4', ('time',)); v[:] = job['ka'](nt); v.units = 'thousand years before present (model time axis)'
        v = d.createVariable('lat', 'f4', ('lat',)); v[:] = lat[i0:i1]; v.units = 'degrees_north'
        v = d.createVariable('lon', 'f4', ('lon',)); v[:] = lon[j0:j1]; v.units = 'degrees_east'
        v = d.createVariable(var, 'f4', ('time', 'lat', 'lon'), zlib=True, complevel=4,
                             chunksizes=(min(CHUNK, nt), i1 - i0, j1 - j0), fill_value=np.float32(-9999.0))
        for a in src[var].ncattrs():
            if a not in ('_FillValue', 'missing_value'):
                try: setattr(v, a, src[var].getncattr(a))
                except Exception: pass
        done = d.createVariable('chunk_done', 'i1', ('time',)); done[:] = 0
        d.source_url = url; d.citation = job['cite']; d.note = job['note']
        d.accessed = f'{datetime.date.today():%Y-%m-%d}'
        d.licence = 'ICCP: "CC BY 4.0 ... CC BY-NC" (contradictory; attribution + noncommercial). Not for redistribution in the GitHub repo.'
        d.close()
    d = Dataset(out, 'a')
    try:
        todo = [s for s in range(0, nt, CHUNK) if not d['chunk_done'][s:min(s + CHUNK, nt)].all()]
        log(f'{name}/{var}: {nt} steps, box lat[{i0}:{i1}] lon[{j0}:{j1}], {len(todo)} chunks to fetch')
        for s in todo:
            e = min(s + CHUNK, nt)
            for attempt in range(6):
                try:
                    a = src[var][s:e, i0:i1, j0:j1]
                    break
                except Exception as ex:
                    log(f'  retry {attempt + 1} {var}[{s}:{e}]: {ex!r}'[:200])
                    time.sleep(10 * (attempt + 1))
                    try: src.close()
                    except Exception: pass
                    src, url = open_remote(job['path'])
            else:
                raise RuntimeError(f'{var}[{s}:{e}] failed after retries')
            a = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(a, dtype=np.float32)), np.float32(-9999.0))
            d[var][s:e] = a; d['chunk_done'][s:e] = 1; d.sync()
            log(f'  {var}[{s}:{e}] ok')
            time.sleep(PAUSE)
    finally:
        d.close(); src.close()
    return out


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def verify(path, var):
    d = Dataset(path)
    ok = bool(d['chunk_done'][:].all())
    a = d[var]
    nt = a.shape[0]
    sample = np.ma.filled(a[nt // 2], np.nan)
    stats = dict(shape=a.shape, complete=ok, mid_min=float(np.nanmin(sample)), mid_max=float(np.nanmax(sample)),
                 ka_first=float(d['ka_bp'][0]), ka_last=float(d['ka_bp'][-1]))
    d.close()
    return stats


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    which = sys.argv[1:] or list(JOBS)
    files = []
    for name in which:
        job = JOBS[name]
        for var in job['vars']:
            p = fetch_var(job, name, var)
            log(f'verify {os.path.basename(p)}: {verify(p, var)}')
            files.append(p)
    allnc = sorted(f for f in os.listdir(OUT) if f.endswith('.nc'))
    with open(os.path.join(OUT, 'SHA256SUMS'), 'w') as f:
        for fn in allnc:
            f.write(f'{sha256(os.path.join(OUT, fn))}  {fn}\n')
    log('done; SHA256SUMS written for ' + ', '.join(allnc))
