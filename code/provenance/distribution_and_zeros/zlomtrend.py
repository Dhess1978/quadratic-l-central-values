"""34.79 odv-lr-zlomtrend podle PRAVIDLA_3479 (doplněk, pravidlo zapsáno před výpočtem). Hybrid jako 34.77 (gen_vddata.py): aritmetika × USp(2N), N z prvního momentu."""
import json, math, sys, numpy as np, mpmath as mp
sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
LH = json.load(open('/home/claude/br76/data/lh_ref.json')); K1 = LH['k1']; a1 = LH['a1']; kap_ar = LH['kap_ar']
lnE = lambda N: 2 * N * math.log(2) + math.lgamma(N + 1.5) - math.lgamma(1.5) - math.lgamma(2 * N + 2) + math.lgamma(N + 2)
def Nk1(mq):
    t = math.log(mq / a1); lo, hi = 0.01, 60.0
    for _ in range(200):
        m = (lo + hi) / 2; lo, hi = (m, hi) if lnE(m) < t else (lo, m)
    return (lo + hi) / 2
ps = cfkrs.primes_upto(2000000).astype(float); xx = ps ** -0.5
ln_a = lambda s: float(np.sum(s * (s + 1) / 2 * np.log1p(-1 / ps) + np.log(1 / (ps + 1) + ps / (ps + 1) * 0.5 * ((1 - xx) ** -s + (1 + xx) ** -s))))
mp.mp.dps = 20; G = lambda z: mp.log(mp.barnesg(z))
lnM = lambda N, s: float(2 * N * s * mp.log(2) + G(2 * N + 2) - G(N + 2) - (G(N + 1.5) - G(1.5)) + (G(N + 1.5 + s) - G(1.5 + s)) - (G(2 * N + 2 + s) - G(N + 2 + s)))
def _Fm(m, z):
    z = mp.mpf(z); return (z - 1) * mp.psi(0, z) - z if m == 1 else (mp.psi(0, z) + (z - 1) * mp.psi(1, z) if m == 2 else 2 * mp.psi(1, z) + (z - 1) * mp.psi(2, z))
def usp_kap(N): return [float((2 * N * mp.log(2) if m == 1 else 0) + _Fm(m, N + 1.5) - _Fm(m, 1.5) - _Fm(m, 2 * N + 2) + _Fm(m, N + 2)) for m in (1, 2, 3)]
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; DV = np.concatenate([v['D'] for v in V]).astype(np.int64); LV = np.concatenate([v['L'] for v in V])
B = np.load('/home/claude/br78/data/D100bc.npz'); DB, LB = B['D'].astype(np.int64), B['L']
PAS = [('V', 1e5, 1e6, 20000), ('V', 1e6, 3e6, 20000), ('V', 3e6, 1e7 + 1, 20000), ('B', 100100000, 102000000, 2000)]
out = {'pasma': []}
for zdroj, lo, hi, w in PAS:
    Dz, Lz = (DV, LV) if zdroj == 'V' else (DB, LB)
    for nm, s, a in (('kl', 1, 0), ('zap', -1, 1)):
        m = (np.sign(Dz) == s) & (np.abs(Dz) >= lo) & (np.abs(Dz) < hi); d = np.abs(Dz[m]); L = Lz[m]; n = len(L)
        N = Nk1(K1['c0'][a] + K1['c1'] * float(np.mean(np.log(d.astype(float))))); b = ((d - int(lo)) // w).astype(np.int64); ub, bi = np.unique(b, return_inverse=True)
        row = dict(zdroj=zdroj, s=nm, lo=lo, hi=hi, n=n, N=N, w=w, bloku=len(ub), rel={}, e={}, dat={}, hyb={})
        for sx in (0.5, -0.5):
            v = L ** sx; mu = float(v.mean()); S = np.bincount(bi, v - mu); e = math.sqrt(float(np.sum(S ** 2))) / (n * mu)
            hyb = math.exp(ln_a(sx) + lnM(N, sx)); row['dat'][str(sx)] = mu; row['hyb'][str(sx)] = hyb; row['rel'][str(sx)] = mu / hyb - 1; row['e'][str(sx)] = e * mu / hyb
        y = np.log(L); c = y - y.mean(); k3 = float(np.mean(c ** 3))
        jk = np.array([float(np.mean((y[bi != t] - y[bi != t].mean()) ** 3)) for t in range(len(ub))]); e3 = math.sqrt((len(jk) - 1) / len(jk) * float(np.sum((jk - jk.mean()) ** 2)))
        k3h = kap_ar[2] + usp_kap(N)[2]; row.update(k3=k3, k3h=k3h, e3=e3, z3=(k3 - k3h) / e3)
        out['pasma'].append(row)
        print('%s %-3s [%.0e, %.0e) n %7d N %.2f  E[L^½]/hyb−1 %+.4f ± %.4f  E[L^−½]/hyb−1 %+.4f ± %.4f  κ₃ %.4f hyb %.4f z %+.1f' % (zdroj, nm, lo, hi, n, N, row['rel']['0.5'], row['e']['0.5'], row['rel']['-0.5'], row['e']['-0.5'], k3, k3h, row['z3']), flush=True)
# verdikt: pásmo [3·10⁶, 10⁷] proti 10⁸
P = {(r['zdroj'], r['s'], r['lo']): r for r in out['pasma']}; zd = []
for nm in ('kl', 'zap'):
    a, b = P[('V', nm, 3e6)], P[('B', nm, 100100000)]
    for sx in ('0.5', '-0.5'):
        d = abs(a['rel'][sx]) - abs(b['rel'][sx]); e = math.hypot(a['e'][sx], b['e'][sx]); zd.append(dict(s=nm, sx=sx, d=d, e=e, z=d / e, rel_mid=a['rel'][sx], rel_top=b['rel'][sx]))
out['srov'] = zd; out['verdikt'] = 'O' if all(x['z'] > 3 for x in zd) else ('Z' if any(x['z'] < -3 for x in zd) else 'N')
# popisně: fit a + b/ln|D| přes 4 pásma (vážený)
fit = {}
for nm in ('kl', 'zap'):
    for sx in ('0.5', '-0.5'):
        R = [r for r in out['pasma'] if r['s'] == nm]; x = np.array([1 / math.log(math.sqrt(r['lo'] * r['hi'])) for r in R]); y = np.array([r['rel'][sx] for r in R]); e = np.array([r['e'][sx] for r in R])
        W = 1 / e ** 2; A = np.vstack([np.ones_like(x), x]).T; C = np.linalg.inv(A.T @ (A * W[:, None])); p = C @ (A.T @ (W * y)); fit['%s_%s' % (nm, sx)] = dict(a=float(p[0]), ea=float(math.sqrt(C[0, 0])), b=float(p[1]), eb=float(math.sqrt(C[1, 1])))
out['fit'] = fit
for x in zd: print('trend', x['s'], x['sx'], 'střed %+.4f → 10⁸ %+.4f  d %+.4f ± %.4f  z %+.1f' % (x['rel_mid'], x['rel_top'], x['d'], x['e'], x['z']))
print('fit', json.dumps(fit)); print('verdikt', out['verdikt'])
json.dump(out, open('/home/claude/br79/data/zlomtrend.json', 'w'), indent=1)
