import sys, json, time, numpy as np
sys.path.insert(0, '/home/claude/br76/py')
import lcheck as lc
rng = np.random.default_rng(20261110)
spf = lc.spf_table(60000)
fam = {F: np.load('/home/claude/br76/data/%s.npz' % F) for F in ['B10', 'C100', 'R10', 'D100', 'V5', 'V10']}
sample = []
for F, k in [('V5', 400), ('V10', 400), ('C100', 100), ('R10', 100), ('D100', 100), ('B10', 100)]:
    D = fam[F]['D']; L = fam[F]['L']; idx = rng.choice(len(D), k, replace=False)
    sample += [(F, int(D[i]), float(L[i])) for i in idx]
allD = np.concatenate([fam[F]['D'] for F in fam]); allL = np.concatenate([fam[F]['L'] for F in fam]); allF = np.concatenate([[F]*len(fam[F]['D']) for F in fam])
for i in np.argsort(allL)[:20]: sample.append(('min:' + allF[i], int(allD[i]), float(allL[i])))
t0 = time.time(); rows = []
for F, D, Lp in sample:
    Lc = lc.L_half(D, spf); rows.append((F, D, Lp, Lc, abs(Lc - Lp), abs(Lc - Lp) / abs(Lp)))
rel = np.array([r[5] for r in rows]); ab = np.array([r[4] for r in rows])
print('n', len(rows), 'time %.0fs' % (time.time() - t0))
print('max rel %.2e  median rel %.2e  max abs %.2e' % (rel.max(), np.median(rel), ab.max()))
for r in sorted(rows, key=lambda r: -r[5])[:6]: print('  worst', r[0], r[1], '%.15g %.15g rel %.1e' % (r[2], r[3], r[5]))
for r in rows[-20:][:8]: print('  small', r[0], r[1], '%.6e %.6e rel %.1e' % (r[2], r[3], r[5]))
json.dump(dict(n=len(rows), max_rel=float(rel.max()), median_rel=float(np.median(rel)), max_abs=float(ab.max()),
               rows=[dict(F=r[0], D=r[1], L_pari=r[2], L_afe=r[3]) for r in rows]), open('/home/claude/br76/data/kontrola_L_scipy.json', 'w'))
