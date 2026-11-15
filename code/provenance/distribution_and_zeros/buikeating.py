"""34.79 odv-lr-split2 a odv-lr-zx podle PRAVIDLA_3479 (pravidla zapsaná před výpočtem).
P_X = exp Σ_{n ≤ X} Λ(n)χ(n)/(√n log n) (Bui–Keating, jen mocniny prvočísel ≤ X), Z_X = L/P_X."""
import json, math, numpy as np, mpmath as mp
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; DV = np.concatenate([v['D'] for v in V]).astype(np.int64); LV = np.concatenate([v['L'] for v in V])
B = np.load('/home/claude/br78/data/D100bc.npz'); DB, LB = B['D'].astype(np.int64), B['L']
def kron(Dv, p):
    if p == 2: r = Dv % 8; return np.where(r % 2 == 0, 0, np.where((r == 1) | (r == 7), 1, -1))
    r = Dv % p; q = np.ones_like(r); e = (p - 1) // 2; b = r.copy()
    while e:
        if e & 1: q = (q * b) % p
        b = (b * b) % p; e >>= 1
    return np.where(r == 0, 0, np.where(q == 1, 1, -1))
PR = [p for p in range(2, 30) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
def lnPX(Dv, X):
    s = np.zeros(len(Dv))
    for p in PR:
        if p > X: break
        c = kron(Dv, p).astype(float); m = 1
        while p ** m <= X: s += c ** m / (m * p ** (m / 2)); m += 1
    return s
mp.mp.dps = 20; G = lambda z: mp.log(mp.barnesg(z))
lnM = lambda N, s: float(2 * N * s * mp.log(2) + G(2 * N + 2) - G(N + 2) - (G(N + 1.5) - G(1.5)) + (G(N + 1.5 + s) - G(1.5 + s)) - (G(2 * N + 2 + s) - G(N + 2 + s)))
EG = math.exp(0.5772156649015329)
PAS = [('V', 1e5, 1e6, 20000), ('V', 3e6, 1e7 + 1, 20000), ('B', 100100000, 102000000, 2000)]
def jk_S(L, PX, ZX, bi, nb, k):
    A = np.bincount(bi, L ** k, nb); Bp = np.bincount(bi, PX ** k, nb); Bz = np.bincount(bi, ZX ** k, nb); n = np.bincount(bi, None, nb)
    f = lambda a, b, c, m: (a / m) / ((b / m) * (c / m)); full = f(A.sum(), Bp.sum(), Bz.sum(), n.sum())
    jk = np.array([f(A.sum() - A[i], Bp.sum() - Bp[i], Bz.sum() - Bz[i], n.sum() - n[i]) for i in range(nb)]); return float(full), float(math.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2)))
out = {'split': [], 'zx': []}
for zdroj, lo, hi, w in PAS:
    Dz, Lz = (DV, LV) if zdroj == 'V' else (DB, LB)
    for nm, sg in (('kl', 1), ('zap', -1)):
        m = (np.sign(Dz) == sg) & (np.abs(Dz) >= lo) & (np.abs(Dz) < hi); d = Dz[m]; L = Lz[m]; ad = np.abs(d)
        b = (ad - int(lo)) // w; ub, bi = np.unique(b, return_inverse=True); nb = len(ub)
        for X in (5, 10, 20):
            PX = np.exp(lnPX(d, X)); ZX = L / PX
            if X == 10:
                for k in (1, 2):
                    S, e = jk_S(L, PX, ZX, bi, nb, k); out['split'].append(dict(zdroj=zdroj, lo=lo, s=nm, k=k, S=S, e=e))
                    print('split %s %-3s [%.0e) k=%d S %.4f ± %.4f' % (zdroj, nm, lo, k, S, e), flush=True)
            if lo >= 3e6:
                Ns = np.log(ad.astype(float)) / (EG * math.log(X)); grid = np.linspace(Ns.min(), Ns.max(), 9)
                for k in (1, 2):
                    mk = np.array([math.exp(lnM(float(N), k)) for N in grid]); pred = float(np.mean(np.interp(Ns, grid, mk)))
                    v = ZX ** k; mu = float(v.mean()); S_ = np.bincount(bi, v - mu, nb); e = math.sqrt(float(np.sum(S_ ** 2))) / len(v)
                    N2 = Ns / 2; grid2 = np.linspace(N2.min(), N2.max(), 9); mk2 = np.array([math.exp(lnM(float(N), k)) for N in grid2]); pred2 = float(np.mean(np.interp(N2, grid2, mk2)))
                    cK = {1: 1 / math.sqrt(2), 2: 1 / 12}[k]; pred_as = float(np.mean(cK * (2 * N2) ** (k * (k + 1) / 2)))
                    out['zx'].append(dict(zdroj=zdroj, lo=lo, s=nm, X=X, k=k, N=float(Ns.mean()), EZ=mu, eEZ=e, usp=pred, pomer=mu / pred, e_pomer=e / pred,
                                          N_bk=float(N2.mean()), usp_bk=pred2, pomer_bk=mu / pred2, asympt_bk=pred_as, pomer_asympt=mu / pred_as))
                    print('zx %s %-3s X=%2d k=%d N %.3f  E[Z^k] %.4f ± %.4f  USp %.4f  poměr %.4f' % (zdroj, nm, X, k, Ns.mean(), mu, e, pred, mu / pred), flush=True)
# verdikty
vs = []
for nm in ('kl', 'zap'):
    for k in (1, 2):
        a = [r for r in out['split'] if r['s'] == nm and r['k'] == k and r['lo'] == 1e5][0]; c = [r for r in out['split'] if r['s'] == nm and r['k'] == k and r['zdroj'] == 'B'][0]
        d = abs(a['S'] - 1) - abs(c['S'] - 1); e = math.hypot(a['e'], c['e']); vs.append(dict(s=nm, k=k, d=d, e=e, z=d / e))
out['split_srov'] = vs; out['split_verdikt'] = 'O' if all(v['z'] > 3 for v in vs) else ('Z' if any(v['z'] < -3 for v in vs) else 'N')
dev = [abs(r['pomer'] - 1) for r in out['zx']]; out['zx_maxdev'] = max(dev); out['zx_verdikt'] = 'O' if max(dev) < 0.05 else ('Z' if max(dev) > 0.15 else 'N')
dev2 = [abs(r['pomer_bk'] - 1) for r in out['zx']]; out['zx_bk_maxdev'] = max(dev2); out['zx_bk_verdikt'] = 'O' if max(dev2) < 0.05 else ('Z' if max(dev2) > 0.15 else 'N')
for r in out['zx']: print('zxBK', r['zdroj'], r['s'], r['X'], r['k'], 'N_bk %.3f poměr přesně %.4f, asympt. %.4f' % (r['N_bk'], r['pomer_bk'], r['pomer_asympt']))
print('zx_bk', round(max(dev2), 4), out['zx_bk_verdikt'])
print('split', out['split_verdikt'], [round(v['z'], 1) for v in vs], '| zx max|poměr−1|', round(max(dev), 4), out['zx_verdikt'])
json.dump(out, open('/home/claude/br79/data/buikeating.json', 'w'), indent=1)
