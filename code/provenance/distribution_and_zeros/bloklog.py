"""34.79 odv-lr-bloklog podle PRAVIDLA_3479 (doplněk; pravidlo zapsáno před výpočtem)."""
import json, math, numpy as np
Q = json.load(open('/home/claude/br76/data/cfkrs_Q_final_v2.json'))
def model(k, a, x):
    c = np.array(Q['k%d_a%d' % (k, a)]['coeffs'], float); return np.polyval(c[::-1], x)
B = np.load('/home/claude/br78/data/D100bc.npz'); DB, LB = B['D'].astype(np.int64), B['L']
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; DV = np.concatenate([v['D'] for v in V]).astype(np.int64); LV = np.concatenate([v['L'] for v in V])
def rho(Dv, Lv, s, a, k, lo, hi, w):
    m = (np.sign(Dv) == s) & (np.abs(Dv) >= lo) & (np.abs(Dv) < hi); d = np.abs(Dv[m]); r = Lv[m] ** k - model(k, a, np.log(d.astype(float)))
    sig = r.std(); b = (d - lo) // w; nb = int(b.max()) + 1; S = np.bincount(b, r, nb); n = np.bincount(b, None, nb); full = n > 0.5 * n.mean()
    z = S[full] / (sig * np.sqrt(n[full])); rr = float(z.std(ddof=1)); return rr, rr / math.sqrt(2 * (int(full.sum()) - 1)), int(full.sum())
LA, LBB = 0.5 * math.log(2e6), 0.5 * math.log(1.01e8); U = [0.6, 0.7, 0.8, 0.9]; out = {'LA': LA, 'LB': LBB, 'U': U, 'srov': []}; zs = []
for nm, s, a in (('kl', 1, 0), ('zap', -1, 1)):
    for k in (1, 2):
        for u in U:
            wA, wB = int(round(math.exp(u * LA))), int(round(math.exp(u * LBB)))
            A = rho(DV, LV, s, a, k, 1000000, 3000000, wA); Bq = rho(DB, LB, s, a, k, 100100000, 102000000, wB)
            z = (A[0] - Bq[0]) / math.hypot(A[1], Bq[1]); zs.append(z)
            out['srov'].append(dict(s=nm, k=k, u=u, wA=wA, wB=wB, rhoA=A[0], errA=A[1], nA=A[2], rhoB=Bq[0], errB=Bq[1], nB=Bq[2], z=z))
            print('%-3s k=%d u=%.1f wA %5d wB %5d  ρA %.4f±%.4f  ρB %.4f±%.4f  z %+.1f' % (nm, k, u, wA, wB, A[0], A[1], Bq[0], Bq[1], z), flush=True)
zm = max(abs(z) for z in zs); out['zmax'] = zm; out['verdikt'] = 'O' if zm < 3 else ('Z' if zm >= 4 else 'N')
# popisně: u_B, při kterém ρ_B(u_B) = ρ_A(u) (interpolace v ln w) — jemná mřížka w_B
pop = {}
for nm, s, a in (('kl', 1, 0), ('zap', -1, 1)):
    for k in (1, 2):
        wg = np.unique(np.round(np.exp(np.linspace(math.log(150), math.log(20000), 25))).astype(int)); rb = np.array([rho(DB, LB, s, a, k, 100100000, 102000000, int(w))[0] for w in wg])
        r = []
        for u in U:
            ra = [c for c in out['srov'] if c['s'] == nm and c['k'] == k and c['u'] == u][0]['rhoA']
            if rb.min() <= ra <= rb.max(): lw = float(np.interp(-ra, -rb, np.log(wg))); r.append(dict(u=u, uB=lw / LBB))
            else: r.append(dict(u=u, uB=None))
        pop['%s_k%d' % (nm, k)] = r
out['popis_uB'] = pop; print('max|z| %.2f verdikt %s' % (zm, out['verdikt'])); print(json.dumps(pop))
json.dump(out, open('/home/claude/br79/data/bloklog.json', 'w'), indent=1)
