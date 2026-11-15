"""Distribution of ln L(1/2, chi_D): data cumulants per |D| range and sign vs hybrid prediction
ln L = ln(arithmetic) + ln Lambda_{USp(2N)}(1), 2N = ln(|D|/pi) (Keating-Snaith identification)."""
import json, math, numpy as np, mpmath as mp, sys
sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]
D = np.concatenate([v['D'] for v in V]); L = np.concatenate([v['L'] for v in V])
# arithmetic cumulants: ln a(s), a(s) = prod_p (1-1/p)^{s(s+1)/2} [1/(p+1) + p/(p+1) ((1-p^-1/2)^-s + (1+p^-1/2)^-s)/2]
ps = cfkrs.primes_upto(2000000).astype(float)
def ln_a(s):
    x = ps ** -0.5
    loc = s * (s + 1) / 2 * np.log1p(-1 / ps) + np.log(1 / (ps + 1) + ps / (ps + 1) * 0.5 * ((1 - x) ** -s + (1 + x) ** -s))
    return float(np.sum(loc))
sg = np.linspace(-0.6, 0.6, 61); va = np.array([ln_a(s) for s in sg])
cf = np.polynomial.polynomial.polyfit(sg, va, 14)
kap_ar = [cf[m] * math.factorial(m) for m in range(1, 5)]
print('arithmetic cumulants k1..k4', kap_ar, ' check ln a(1)=%.10f vs a_1 %.10f' % (ln_a(1.0), math.log(cfkrs.a_k(1, 2000000))))
mp.mp.dps = 30
def usp_cum(N):
    N = mp.mpf(N)
    def lnM(s):
        return 2 * N * s * mp.log(2) + mp.log(mp.barnesg(2 * N + 2)) - mp.log(mp.barnesg(N + 2)) - (mp.log(mp.barnesg(N + mp.mpf(3) / 2)) - mp.log(mp.barnesg(mp.mpf(3) / 2))) \
            + (mp.log(mp.barnesg(N + mp.mpf(3) / 2 + s)) - mp.log(mp.barnesg(mp.mpf(3) / 2 + s))) - (mp.log(mp.barnesg(2 * N + 2 + s)) - mp.log(mp.barnesg(N + 2 + s)))
    return [float(mp.diff(lnM, 0, m)) for m in range(1, 5)]
# sanity: integer N, direct product formula for E[Lambda^1] and E[Lambda^2]
for N in (3, 6):
    direct = lambda s: 2 ** (2 * N * s) * mp.fprod([mp.gamma(1 + N + j) * mp.gamma(mp.mpf(1) / 2 + s + j) / (mp.gamma(mp.mpf(1) / 2 + j) * mp.gamma(1 + s + N + j)) for j in range(1, N + 1)])
    print('USp check N=%d' % N, float(mp.log(direct(1))), float(mp.diff(lambda s: mp.log(direct(s)), 0, 1)), usp_cum(N)[0])
groups = [(1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 3e6), (3e6, 1e7 + 1)]
out = []
for s, name in ((1, 'kl'), (-1, 'zap')):
    for lo, hi in groups:
        m = (np.sign(D) == s) & (np.abs(D) >= lo) & (np.abs(D) < hi)
        y = np.log(L[m]); n = len(y)
        mu = y.mean(); c = y - mu; var = (c ** 2).mean(); k3 = (c ** 3).mean(); k4 = (c ** 4).mean() - 3 * var ** 2
        med = float(np.median(np.abs(D[m]))); N = math.log(med / math.pi) / 2
        u = usp_cum(N)
        pred = [kap_ar[i] + u[i] for i in range(4)]
        out.append(dict(sign=name, lo=lo, hi=hi, n=n, absD_med=med, N=N, data=[mu, var, k3, k4], usp=u, arith=kap_ar, pred=pred,
                        skew=k3 / var ** 1.5, exkurt=k4 / var ** 2, pred_skew=pred[2] / pred[1] ** 1.5, pred_exkurt=pred[3] / pred[1] ** 2))
        print('%s |D| %.0e-%.0e n=%7d N=%.2f  mean %.4f/%.4f  var %.4f/%.4f  k3 %.4f/%.4f  k4 %.4f/%.4f' % (name, lo, hi, n, N, mu, pred[0], var, pred[1], k3, pred[2], k4, pred[3]), flush=True)
json.dump(dict(groups=out, kap_ar=kap_ar), open('/home/claude/br76/data/rozdeleni.json', 'w'), indent=1)
# ---- refinement: N per sign fixed by the first moment, a_1 * E_USp(2N)[Lambda] = Q_1(ln|D|) (no free parameter)
Q1 = {0: (-0.4889851880, 0.3522211004), 1: (0.6175500334, 0.3522211004)}   # const, slope (P=1e7 analytic)
a1 = math.exp(ln_a(1.0))
def lnM1(N):
    N = mp.mpf(N); s = 1
    return float(2 * N * s * mp.log(2) + mp.log(mp.barnesg(2 * N + 2)) - mp.log(mp.barnesg(N + 2)) - (mp.log(mp.barnesg(N + mp.mpf(3) / 2)) - mp.log(mp.barnesg(mp.mpf(3) / 2)))
                 + (mp.log(mp.barnesg(N + mp.mpf(3) / 2 + s)) - mp.log(mp.barnesg(mp.mpf(3) / 2 + s))) - (mp.log(mp.barnesg(2 * N + 2 + s)) - mp.log(mp.barnesg(N + 2 + s))))
out2 = []
for s, a, name in ((1, 0, 'kl'), (-1, 1, 'zap')):
    for lo, hi in groups:
        m = (np.sign(D) == s) & (np.abs(D) >= lo) & (np.abs(D) < hi)
        y = np.log(L[m]); mu = y.mean(); c = y - mu; var = (c ** 2).mean(); k3 = (c ** 3).mean(); k4 = (c ** 4).mean() - 3 * var ** 2
        x = np.log(np.abs(D[m]).astype(float)); target = math.log(np.mean(Q1[a][0] + Q1[a][1] * x) / a1)
        N = float(mp.findroot(lambda NN: lnM1(NN) - target, 5.0))
        u = usp_cum(N); pred = [kap_ar[i] + u[i] for i in range(4)]
        # bootstrap-free rough errors: blocks of 1000 consecutive |D|
        out2.append(dict(sign=name, lo=lo, hi=hi, n=int(m.sum()), N_from_k1=N, data=[mu, var, k3, k4], pred=pred))
        print('k1-matched %s |D| %.0e-%.0e N=%.2f  mean %.4f/%.4f  var %.4f/%.4f  k3 %.4f/%.4f  k4 %.4f/%.4f' % (name, lo, hi, N, mu, pred[0], var, pred[1], k3, pred[2], k4, pred[3]), flush=True)
json.dump(dict(groups=out, groups_k1=out2, kap_ar=kap_ar), open('/home/claude/br76/data/rozdeleni.json', 'w'), indent=1)
