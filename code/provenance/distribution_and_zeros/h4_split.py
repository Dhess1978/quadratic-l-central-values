"""34.78 H4 odv-lq-split: E[L^k] = E[P_X^k]·E[Z_X^k], P_X = Π_{p≤X}(1 − χ_D(p)/√p)⁻¹, Z_X = L/P_X; |D| ∈ [3·10⁶, 10⁷]."""
import json, math, numpy as np
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; D = np.concatenate([v['D'] for v in V]).astype(np.int64); L = np.concatenate([v['L'] for v in V])
def chi(Dv, p):
    if p == 2: r = Dv % 8; return np.where(Dv % 2 == 0, 0, np.where((r == 1) | (r == 7), 1, -1))
    r = Dv % p; q = np.ones_like(r)
    e = (p - 1) // 2; b = r.copy()
    while e:            # rychlé umocňování r^((p−1)/2) mod p
        if e & 1: q = (q * b) % p
        b = (b * b) % p; e >>= 1
    return np.where(r == 0, 0, np.where(q == 1, 1, -1))
P = [p for p in range(2, 301) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
out = {}
for nm, s in (('kl', 1), ('zap', -1)):
    m = (np.sign(D) == s) & (np.abs(D) >= 3e6); d = D[m]; l = L[m]; blk = (np.abs(d) - 1) // 20000; ub, bi = np.unique(blk, return_inverse=True); nbk = len(ub)
    lnP = np.zeros(len(d)); o = {}
    for p in P:
        lnP -= np.log1p(-chi(d, p) / math.sqrt(p))
        if p in (7, 29, 97, 293):   # X = 10, 30, 100, 300 (poslední prvočíslo ≤ X)
            X = {7: 10, 29: 30, 97: 100, 293: 300}[p]; Pp = np.exp(lnP); Z = l / Pp; r = {}
            for k in (1, 2, 3):
                A = np.bincount(bi, l ** k, nbk); Bp = np.bincount(bi, Pp ** k, nbk); Bz = np.bincount(bi, Z ** k, nbk); n = np.bincount(bi, None, nbk)
                def S(a, b, c, nn): return (a.sum() / nn.sum()) / ((b.sum() / nn.sum()) * (c.sum() / nn.sum()))
                full = S(A, Bp, Bz, n); jk = np.array([S(A.sum() - A[i] + 0 * A, Bp.sum() - Bp[i] + 0 * Bp, Bz.sum() - Bz[i] + 0 * Bz, n.sum() - n[i] + 0 * n) if False else ((A.sum() - A[i]) / (n.sum() - n[i])) / (((Bp.sum() - Bp[i]) / (n.sum() - n[i])) * ((Bz.sum() - Bz[i]) / (n.sum() - n[i]))) for i in range(nbk)])
                e = math.sqrt((nbk - 1) / nbk * np.sum((jk - jk.mean()) ** 2)); r[k] = dict(S=float(full), e=e)
            o[X] = dict(S=r, corr=float(np.corrcoef(np.log(Pp), np.log(Z))[0, 1]))
            print(nm, 'X', X, ' '.join('S%d %.4f±%.4f' % (k, r[k]['S'], r[k]['e']) for k in (1, 2, 3)), ' corr(ln P, ln Z) %.3f' % o[X]['corr'], flush=True)
    out[nm] = o
dev = [abs(out[nm][X]['S'][k]['S'] - 1) for nm in out for X in (10, 30, 100) for k in (1, 2, 3)]
out['verdikt'] = 'O' if max(dev) < 0.01 else ('Z' if max(dev) > 0.03 else 'N'); out['max_odchylka'] = max(dev); print('max |S−1| (X ≤ 100)', round(max(dev), 4), 'verdikt', out['verdikt'])
json.dump(out, open('/home/claude/br78/data/h4_split.json', 'w'), indent=1)
