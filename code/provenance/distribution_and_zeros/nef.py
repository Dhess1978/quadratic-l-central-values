"""34.79 odv-lr-nef podle PRAVIDLA_3479 (doplněk, pravidlo zapsáno před výpočtem).
N_ef = argmin Σ (W_USp(2N) − W_poměry)² na koších 0,05 v x ∈ [0, 1]; domněnka poměrů jako 34.62 (kx_pomery.py: Conrey–Snaith s A_D, rozvinutí θ_D);
předpověď střední x₁ a Δ₁ (vážení Λ(1)) z USp(2N_ef) přesně (Jacobiho soubor, tabulka 34.78 h2_vaha.json, lineárně mezi celými N)."""
import sys, os, json, math, csv, time, numpy as np, mpmath as mp, scipy.special as sp
from scipy.optimize import minimize_scalar
t0 = time.time()
DT = 0.002; TG = np.round(np.arange(0, 8.0 + 1e-9, DT), 4); NT = len(TG)
Mp = 200000; s = np.ones(Mp + 1, bool); s[:2] = False
for p in range(2, int(Mp ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0].astype(float); lp = np.log(P)
A1 = np.empty(NT, complex); A2 = np.empty(NT, complex)
for i, tt in enumerate(TG):
    z = np.exp(2j * tt * lp); A1[i] = np.exp(np.sum(np.log(1 - (z - 1) / (P * P - 1)))); A2[i] = np.sum(lp / ((P + 1) * (P * z - 1)))
mp.mp.dps = 30; ZL = np.zeros(NT, complex); ZR = np.zeros(NT, complex)
for i, tt in enumerate(TG):
    if i == 0: continue
    s1 = mp.mpc(1, 2 * float(tt)); ZL[i] = complex(mp.zeta(s1, derivative=1) / mp.zeta(s1)); ZR[i] = complex(mp.zeta(mp.mpc(1, -2 * float(tt))))
GR = {k: np.exp(sp.loggamma(0.25 + 0.5 * k - 0.5j * TG) - sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
IMLG = {k: np.imag(sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
Z1 = ZL + A2; Z2 = {k: GR[k] * ZR * A1 for k in (0, 1)}
print('mřížka', round(time.time() - t0), 's', flush=True)
XB, NB = 0.05, 20
def predl(l, k):
    c = np.real(Z1 - np.exp(-1j * TG * l) * Z2[k]) / math.pi; c[0] = 3 * c[1] - 3 * c[2] + c[3]
    C = np.concatenate([[0], np.cumsum(0.5 * (c[1:] + c[:-1]) * DT)]); th = (0.5 * TG * l + IMLG[k]) / math.pi; te = np.interp(np.arange(NB + 1) * XB, th, TG)
    return 1 + np.diff(np.interp(te, TG, C)) / XB
xs = np.linspace(0, NB * XB, NB * 200 + 1); xs[0] = 1e-9
def W_N(x, N): y = np.pi * x / N; return 1 + 1 / (2 * N) - np.sin((2 * N + 1) * y) / (2 * N * np.sin(np.minimum(y, np.pi - 1e-12)))
def binavg_N(N): v = W_N(xs, N); seg = 0.5 * (v[1:] + v[:-1]) * np.diff(xs); return seg.reshape(NB, 200).sum(1) / XB
H2 = json.load(open('/home/claude/br78/data/h2_vaha.json')); U = {int(N): v for N, v in H2['usp'].items()}
def usp_int(N, j):
    i = min(max(int(math.floor(N)), 1), max(U) - 1); t = N - i; return U[i][j] * (1 - t) + U[i + 1][j] * t
def lngamma_c(z): return sp.loggamma(z)
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as fh:
    rd = csv.reader(fh); next(rd)
    for F, d0, g1, g2, nn in rd: rows.setdefault(F, []).append((int(d0), float(g1)))
out = {'skupiny': {}}; zs = []
for F in ('B10', 'C100', 'R10', 'D100'):
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        pr = [(d0, g1) for d0, g1 in rows[F] if (d0 > 0) == (sg > 0) and d0 in lv and math.isfinite(g1)]
        Dv = np.array([p[0] for p in pr], np.int64); t = np.array([p[1] for p in pr]); x = ((t / 2) * np.log(np.abs(Dv) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi
        n = len(x); mx1 = float(x.mean()); se = float(x.std(ddof=1) / math.sqrt(n))
        l = np.log(np.abs(Dv) / math.pi); edges = np.linspace(l.min(), l.max() + 1e-9, 41); hcnt, _ = np.histogram(l, edges); cen = [l[(l >= edges[i]) & (l < edges[i + 1])].mean() for i in range(40) if hcnt[i] > 0]
        wts = hcnt[hcnt > 0] / n; Wp = sum(w * predl(c, a) for w, c in zip(wts, cen))
        Wp1 = predl(float(l.mean()), a)   # kontrola: předpověď ve středu l
        res = minimize_scalar(lambda N: float(np.sum((binavg_N(N) - Wp) ** 2)), bounds=(1.0, 11.5), method='bounded', options=dict(xatol=1e-6)); Nef = float(res.x)
        mx_u = usp_int(Nef, 0); d_u = usp_int(Nef, 1) - mx_u
        g = H2['skupiny']['%s_%s' % (F, nm)]; zm = (mx1 - mx_u) / se; zd = (g['d'][0] - d_u) / g['e'][0]; zs += [zm, zd]
        out['skupiny']['%s_%s' % (F, nm)] = dict(n=n, l_mean=float(l.mean()), Nef=Nef, res=float(res.fun), mx1=mx1, se=se, mx1_usp=mx_u, z_mx1=zm, d1=g['d'][0], e1=g['e'][0], d1_usp=d_u, z_d1=zd,
                                                 N_x1=g['N_x1'], W_pom=Wp.tolist(), W_usp=binavg_N(Nef).tolist(), max_rozdil_l=float(np.max(np.abs(Wp - Wp1))))
        print('%-9s n %6d N_ef %.3f (N_x1 %.2f)  x̄₁ %.4f ± %.4f USp %.4f z %+.1f | Δ₁ %.4f ± %.4f USp %.4f z %+.1f | rezid %.2e' % ('%s_%s' % (F, nm), n, Nef, g['N_x1'], mx1, se, mx_u, zm, g['d'][0], g['e'][0], d_u, zd, res.fun), flush=True)
zmax = max(abs(z) for z in zs); out['zmax'] = zmax; out['verdikt'] = 'O' if zmax < 3 else ('Z' if zmax >= 4 else 'N')
print('max|z| %.1f verdikt %s' % (zmax, out['verdikt']), round(time.time() - t0), 's')
json.dump(out, open('/home/claude/br79/data/nef.json', 'w'), indent=1)
