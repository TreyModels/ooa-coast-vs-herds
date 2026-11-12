"""Build the 20-arc-minute study grid (lon 25-160E, lat 20S-45N) from ETOPO 10'."""
import numpy as np
from scipy import ndimage

g = np.load('/home/claude/work/etopo10.npy')            # rows: lat 90 -> -90 step 1/6
lat10 = 90 - np.arange(g.shape[0]) / 6
lon10 = -180 + np.arange(g.shape[1]) / 6

LON0, LON1, LAT0, LAT1 = 25, 160, -20, 45
r = np.where((lat10 <= LAT1) & (lat10 > LAT0))[0]
c = np.where((lon10 >= LON0) & (lon10 < LON1))[0]
sub = g[r[0]:r[-1] + 1, c[0]:c[-1] + 1].astype(np.float64)
sub = sub[: sub.shape[0] // 2 * 2, : sub.shape[1] // 2 * 2]
ny, nx = sub.shape[0] // 2, sub.shape[1] // 2
blk = sub.reshape(ny, 2, nx, 2)
# keep the HIGHEST of the 4 sub-cells so narrow land bridges / islands aren't lost,
# and the median for "typical" elevation
zmax = blk.max(axis=(1, 3))
zmed = np.median(blk.reshape(ny, 2, nx, 2).transpose(0, 2, 1, 3).reshape(ny, nx, 4), axis=2)
lat = lat10[r[0]] - (np.arange(ny) * 2 + 0.5) / 6
lon = lon10[c[0]] + (np.arange(nx) * 2 + 0.5) / 6
zmin = blk.min(axis=(1, 3))
np.savez('/home/claude/work/grid.npz', zmin=zmin, zmax=zmax, zmed=zmed, lat=lat, lon=lon)
print(ny, nx, lat[[0, -1]], lon[[0, -1]])
