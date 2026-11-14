import numpy as np, pickle, json, itertools
R = pickle.load(open('wave2.pkl','rb')); A1 = np.load('arr_L.npy')   # wave-1 arrival BP: [run, strat, cell]
S = ['coastal','terrestrial','mixed']; P = list(R[0]['mixed'].keys())
out = {}
# ---------- (1) split times between population pairs
pairs = [('New Guinea (Papuans)','North China (Tianyuan)'),('New Guinea (Papuans)','South India (AASI)'),
         ('North China (Tianyuan)','South India (AASI)'),('New Guinea (Papuans)','Andamans (Onge)'),
         ('New Guinea (Papuans)','N Australia (Aboriginal)'),('Levant (West Eurasian side)','New Guinea (Papuans)')]
def split(a, b, t0):
    if a is None or b is None: return None
    pa, pb = a['path'], b['path']; n = 0
    while n < min(len(pa), len(pb)) and pa[n] == pb[n]: n += 1
    return float(a['t'][n-1]) if n else t0
east = ['South India (AASI)','Andamans (Onge)','Malay Peninsula (Semang)','North China (Tianyuan)','New Guinea (Papuans)']
out['splits'] = {}; out['star'] = {}; out['topology'] = {}
for s in S:
    sp = {f'{a} | {b}': [] for a, b in pairs}; star = []; topo = {'(PNG,China),India': 0, '(PNG,India),China': 0, '(China,India),PNG': 0, 'tie': 0}
    for r in R:
        d = r[s]
        for a, b in pairs:
            v = split(d[a], d[b], r['t0'])
            if v is not None: sp[f'{a} | {b}'].append(v)
        ts = [split(d[a], d[b], r['t0']) for a, b in itertools.combinations(east, 2)]
        if all(t is not None for t in ts): star.append((max(ts) - min(ts)))
        pc, pi, ci = (split(d['New Guinea (Papuans)'], d['North China (Tianyuan)'], r['t0']),
                      split(d['New Guinea (Papuans)'], d['South India (AASI)'], r['t0']),
                      split(d['North China (Tianyuan)'], d['South India (AASI)'], r['t0']))
        if None in (pc, pi, ci): continue
        m = max(pc, pi, ci)   # most recent split = largest BP? no: most recent = smallest BP
        mn = min(pc, pi, ci)
        if [pc, pi, ci].count(mn) > 1: topo['tie'] += 1
        elif mn == pc: topo['(PNG,China),India'] += 1
        elif mn == pi: topo['(PNG,India),China'] += 1
        else: topo['(China,India),PNG'] += 1
    out['splits'][s] = {k: [round(float(np.percentile(v, q))/1000, 1) for q in (10, 50, 90)] for k, v in sp.items() if v}
    out['star'][s] = [round(float(np.percentile(star, q))/1000, 1) for q in (10, 50, 90)]
    tot = sum(topo.values()); out['topology'][s] = {k: round(v / tot, 2) for k, v in topo.items()}
# ---------- (2) serial founder drift: sum of 1/front-width along the founding path
out['drift'] = {}
for s in S:
    F = {p: [] for p in P}
    for r in R:
        for p in P:
            d = r[s][p]
            if d is not None: F[p].append(float(np.sum(1.0 / d['w'])))
    out['drift'][s] = {p: round(float(np.median(v)), 1) for p, v in F.items() if v}
# relative heterozygosity H/H_hub = exp(-c F), c set so mixed-strategy Papuans lose 15% (about the observed
# gap between Papuans and West/South Asians in genome-wide heterozygosity)
c = -np.log(0.85) / out['drift']['mixed']['New Guinea (Papuans)']
out['H_rel'] = {s: {p: round(float(np.exp(-c * f)), 3) for p, f in out['drift'][s].items()} for s in S}
# ---------- (3) ghost ancestry from an earlier wave
ms = [1e-4, 3e-4, 1e-3, 3e-3]
targets = ['New Guinea (Papuans)','N Australia (Aboriginal)','Luzon (Ayta)','Andamans (Onge)','South India (AASI)','North China (Tianyuan)']
out['ghost'] = {}; out['ghost_cells'] = {}
for k1, s1 in enumerate(S):            # wave-1 strategy
    for s2 in S:                       # wave-2 strategy
        key = f'{s1}>{s2}'; G = {t: {str(m): [] for m in ms} for t in targets}; cells = {t: [] for t in targets}
        for i, r in enumerate(R):
            a1 = A1[i % len(A1), k1]
            for t in targets:
                d = r[s2][t]
                if d is None: continue
                occ = a1[d['path']] > d['t']        # wave 1 got there first (and before the snapshot)
                n = int(occ.sum()); cells[t].append(n)
                for m in ms: G[t][str(m)].append(1 - (1 - m) ** n)
        out['ghost'][key] = {t: {mm: round(float(np.median(v)) * 100, 2) for mm, v in G[t].items()} for t in targets}
        out['ghost_cells'][key] = {t: int(np.median(v)) for t, v in cells.items()}
# m that keeps Papuan ghost ancestry at 2% given the median number of wave-1 cells crossed
out['m_star'] = {k: (round(float(1 - 0.98 ** (1 / v['New Guinea (Papuans)'])), 6) if v['New Guinea (Papuans)'] > 0 else None) for k, v in out['ghost_cells'].items()}
json.dump(out, open('gen.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ['topology','star','drift','H_rel','ghost_cells','m_star']}, indent=1))
print(json.dumps(out['ghost']['coastal>coastal'], indent=0)[:600]); print(json.dumps(out['splits']['mixed'], indent=0))
