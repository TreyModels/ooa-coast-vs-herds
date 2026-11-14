import numpy as np, model as M, erectus as E, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, LinearSegmentedColormap
D = np.load('erectus_reach.npz'); n = int(D['n'])
canoe = D['mixed|canoe'] / n; raft = D['mixed|raft'] / n; walk = D['mixed|walk'] / n
r = (M.LAT <= 45) & (M.LAT >= -20); c = (M.LON >= 30) & (M.LON <= 160)
ext = [M.LON[c][0], M.LON[c][-1], M.LAT[r][-1], M.LAT[r][0]]
land_lo = ~M.water_mask(-80.0); land_now = ~M.water_mask(0.0)
rgb = np.ones(land_now.shape + (3,)) * np.array(to_rgb('#e9eef3'))
rgb[land_lo] = to_rgb('#d8d3c8'); rgb[land_now] = to_rgb('#ece8e0')
fig = plt.figure(figsize=(9, 4.5)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
ax.imshow(rgb[np.ix_(r, c)], extent=ext, interpolation='nearest')
def layer(prob, col, a=0.85):
    p = prob.reshape(M.NY, M.NX)[np.ix_(r, c)]
    im = np.zeros(p.shape + (4,)); im[..., :3] = to_rgb(col); im[..., 3] = np.where(p >= 0.2, a * np.clip(p, 0.35, 1), 0)
    ax.imshow(im, extent=ext, interpolation='nearest')
only_canoe = np.clip(canoe - raft, 0, 1)         # places canoes add beyond rare rafting
only_raft = np.clip(raft - walk, 0, 1)           # places rafting adds beyond walking
layer(only_raft, '#6250d6'); layer(only_canoe, '#e34948')
for name, (la, lo, age, rr) in E.PRESENT.items():
    ax.plot(lo, la, 'o', ms=5, mfc='#0b0b0b', mec='white', mew=1)
    ax.annotate(name.split(',')[0], (lo, la), xytext=(4, 3), textcoords='offset points', fontsize=7)
for name, (la, lo) in {'Sahul': (-6, 140), 'Andamans': (12, 92.8), 'Socotra': (12.5, 53.9)}.items():
    ax.plot(lo, la, 'x', ms=7, mew=2, color='#0b0b0b'); ax.annotate(name + ' (none found)', (lo, la), xytext=(5, -9), textcoords='offset points', fontsize=7)
fig.savefig('map_erectus.png', dpi=100); plt.close(fig)
print('canoe-only share in Sahul cells', round(float(canoe[E.ABSENT['Sahul (New Guinea + Australia)']].mean()),2))
