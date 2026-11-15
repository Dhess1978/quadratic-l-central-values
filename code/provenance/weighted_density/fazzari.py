"""34.79 odv-lr-vaha podle PRAVIDLA_3479: vážená (L(½)) hustota nízkých nul z domněnky poměrů (Fazzari 2109.07244, rov. 4.8, K = 2, Q = 1, A → 0)
proti datům studie III (rozdíl vážené − nevážené hustoty, koše 0,2 v x ∈ [0, 1], chyby z h3_hust.json = Poissonův bootstrap 300×).
Na každé D: S(A, D) = Q1 + (d/π)^{−A} g(½+A) Q2 + (d/π)^{−D} g(½+D) Q3 + (d/π)^{−A−D} g(½+A+D) Q4, A → 0 průměrem po kružnici |A| = h (S je v A analytická);
M = ⟨L(½)⟩ = ζ(1+2A) A₁(A) + (d/π)^{−A} g(½+A) ζ(1−2A) A₁(−A) (A → 0). Hustota v t: θ′/π·M + Re S(it)/π. g_S pro sudá D (a = 0), pro lichá Γ-faktor s a = 1.
Nevážená hustota: domněnka poměrů jako 34.62 (kx_pomery.py). Pravidla zapsaná před výpočtem."""
import sys, json, math, csv, time, numpy as np, mpmath as mp, scipy.special as sp
from scipy import stats
t0 = time.time(); mp.mp.dps = 30
Mp = 200000; s = np.ones(Mp + 1, bool); s[:2] = False
for p in range(2, int(Mp ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0].astype(float); lp = np.log(P); RP = P / (P + 1)
pw = lambda e: np.exp(-e * lp)            # p^{−e}, e komplexní skalár
def A1(A):
    x = pw(0.5 + A); E0 = 0.5 * (1 / (1 - x) + 1 / (1 + x)) - 1
    return np.exp(np.sum(np.log(1 - pw(1 + 2 * A)) + np.log(1 + RP * E0)))
def A2(A, C, D):
    x, y, w = pw(0.5 + A), pw(0.5 + C), pw(0.5 + D); a2, c2, ac, ad, cd = pw(1 + 2 * A), pw(1 + 2 * C), pw(1 + A + C), pw(1 + A + D), pw(1 + C + D)
    f1 = np.log(1 - a2) + np.log(1 - c2) + np.log(1 - ac) - np.log(1 - ad) - np.log(1 - cd)
    pp, mm = (1 - x) * (1 - y), (1 + x) * (1 + y); E0 = 0.5 * (1 / pp + 1 / mm) - 1; O = 0.5 * (1 / pp - 1 / mm); br = 1 + RP * (E0 - w * O)
    val = np.exp(np.sum(f1 + np.log(br)))
    d1 = 2 * lp * c2 / (1 - c2) + lp * ac / (1 - ac) - lp * cd / (1 - cd)
    dE0 = 0.5 * (1 / ((1 - x) * (1 - y) ** 2) - 1 / ((1 + x) * (1 + y) ** 2)); dO = 0.5 * (1 / ((1 - x) * (1 - y) ** 2) + 1 / ((1 + x) * (1 + y) ** 2))
    dbr = RP * (dE0 - w * dO) * (-lp * y)
    return val, val * np.sum(d1 + dbr / br)
Z = lambda z: complex(mp.zeta(mp.mpc(z.real, z.imag))); Zd = lambda z: complex(mp.zeta(mp.mpc(z.real, z.imag), derivative=1))
def gS(z, a): return complex(mp.gamma(mp.mpc((1 - z.real + a) / 2, -z.imag / 2)) / mp.gamma(mp.mpc((z.real + a) / 2, z.imag / 2)))   # a = 0: Γ((1−z)/2)/Γ(z/2); a = 1: Γ((2−z)/2)/Γ((1+z)/2)
NC, H = 16, 0.0037; AJ = [H * complex(math.cos(2 * math.pi * (j + 0.5) / NC), math.sin(2 * math.pi * (j + 0.5) / NC)) for j in range(NC)]
DT = 0.002; TG = np.round(np.arange(DT, 1.40 + 1e-9, DT), 4); NT = len(TG)
PC = {}   # (a) → pole kusů: [NT, NC, 4] a M kusy [NC, 2]
for a in (0, 1):
    Mk = np.array([[Z(1 + 2 * A) * A1(A), gS(0.5 + A, a) * Z(1 - 2 * A) * A1(-A)] for A in AJ])
    Q = np.zeros((NT, NC, 4), complex)
    for i, t in enumerate(TG):
        D = 1j * float(t); zd2 = Z(1 + 2 * D); lz = Zd(1 + 2 * D) / zd2; zm2 = Z(1 - 2 * D)
        for j, A in enumerate(AJ):
            vA, dA = A2(A, D, D); vmA, dmA = A2(-A, D, D); v3, _ = A2(A, -D, D); v4, _ = A2(-A, -D, D)
            za, zma = Z(1 + 2 * A), Z(1 - 2 * A); zpd, zmd = Z(1 + A + D), Z(1 - A + D)
            q1 = vA * za / zpd * (lz * zpd + Zd(1 + A + D)) + dA * za
            q2 = vmA * zma / zmd * (lz * zmd + Zd(1 - A + D)) + dmA * zma
            q3 = -v3 * za * zm2 * Z(1 + A - D) / zpd
            q4 = -v4 * zma * zm2 * Z(1 - A - D) / zmd
            Q[i, j] = [q1, gS(0.5 + A, a) * q2, gS(0.5 + D, a) * q3, gS(0.5 + A + D, a) * q4]
        if i % 100 == 0: print('a', a, 't', t, round(time.time() - t0), 's', flush=True)
    PC[a] = (Q, Mk)
np.savez('/home/claude/br79/data/fazzari_kusy.npz', TG=TG, Q0=PC[0][0], Q1=PC[1][0], M0=PC[0][1], M1=PC[1][1], AJ=np.array(AJ))
AJa = np.array(AJ)
def vahova(l, a):
    """pro jedno l = ln(|d|/π): M a Re S(it)/π na TG"""
    Q, Mk = PC[a]; M = np.mean(Mk[:, 0] + np.exp(-AJa * l) * Mk[:, 1])
    E = np.stack([np.ones((NT, NC)), np.exp(-AJa * l)[None, :].repeat(NT, 0), np.exp(-1j * TG * l)[:, None].repeat(NC, 1), np.exp(-(AJa[None, :] + 1j * TG[:, None]) * l)], -1)
    S = np.mean(np.sum(Q * E, -1), 1); return M, np.real(S) / math.pi
# nevážená (34.62)
A1g = np.empty(NT, complex); A2g = np.empty(NT, complex)
for i, tt in enumerate(TG):
    z = np.exp(2j * tt * lp); A1g[i] = np.exp(np.sum(np.log(1 - (z - 1) / (P * P - 1)))); A2g[i] = np.sum(lp / ((P + 1) * (P * z - 1)))
ZLg = np.array([Zd(complex(1, 2 * t)) / Z(complex(1, 2 * t)) for t in TG]); ZRg = np.array([Z(complex(1, -2 * t)) for t in TG])
GRg = {k: np.exp(sp.loggamma(0.25 + 0.5 * k - 0.5j * TG) - sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
IMLG = {k: np.imag(sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
def nevazena(l, a): return np.real(ZLg + A2g - np.exp(-1j * TG * l) * GRg[a] * ZRg * A1g) / math.pi
def kose(l, a, c, edges):
    """∫ c dt přes koše v x = θ/π (θ = ½ t l + Im lnΓ(¼ + a/2 + it/2)); c na TG, v t = 0 extrapolace"""
    tt = np.concatenate([[0], TG]); c0 = 3 * c[0] - 3 * c[1] + c[2]; cc = np.concatenate([[c0], c]); C = np.concatenate([[0], np.cumsum(0.5 * (cc[1:] + cc[:-1]) * DT)])
    th = (0.5 * tt * l + np.concatenate([[0], IMLG[a]])) / math.pi; te = np.interp(edges, th, tt); assert th[-1] > edges[-1], (l, a, th[-1])
    return np.diff(np.interp(te, tt, C))
EDG = np.linspace(0, 1, 6); BW = 0.2
H3 = json.load(open('/home/claude/br78/data/h3_hust.json'))
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as fh:
    rd = csv.reader(fh); next(rd)
    for F, d0, g1, g2, nn in rd: rows.setdefault(F, []).append(int(d0))
out = {'skupiny': {}, 'h': H, 'nc': NC, 'dt': DT}; ps = []
for F in ('B10', 'C100', 'R10', 'D100'):
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        Dv = np.array([d for d in rows[F] if (d > 0) == (sg > 0)]); l = np.log(np.abs(Dv) / math.pi)
        ed = np.linspace(l.min(), l.max() + 1e-9, 21); cnt, _ = np.histogram(l, ed); cen = [l[(l >= ed[i]) & (l < ed[i + 1])].mean() for i in range(20) if cnt[i] > 0]; wts = cnt[cnt > 0].astype(float)
        num = np.zeros(5); den = 0.0; nu = np.zeros(5); Ms = []
        for w, lc in zip(wts, cen):
            M, cw = vahova(lc, a); Ms.append(M.real); num += w * (M.real * BW + kose(lc, a, cw, EDG)); den += w * M.real * BW
            nu += w * (BW + kose(lc, a, nevazena(lc, a), EDG))
        vaz = num / den; nev = nu / (wts.sum() * BW); pred = vaz - nev
        g = H3['skupiny']['%s_%s' % (F, nm)]; dat = np.array(g['rozdil'][:5]); e = np.array(g['chyba'][:5]); zb = (dat - pred) / e; c2 = float(np.sum(zb ** 2)); p = float(stats.chi2.sf(c2, 5)); ps.append(p)
        usp = np.array(g['usp'][:5])
        out['skupiny']['%s_%s' % (F, nm)] = dict(n=int(len(Dv)), l=float(l.mean()), M=float(np.average(Ms, weights=wts)), pred=pred.tolist(), vazena_pred=vaz.tolist(), nevazena_pred=nev.tolist(),
            nevazena_dat=g['nevazena'][:5], vazena_dat=g['vazena'][:5], rozdil=dat.tolist(), chyba=e.tolist(), z=zb.tolist(), chi2=c2, p=p, usp=usp.tolist(), chi2_usp=float(np.sum(((dat - usp) / e) ** 2)))
        print('%-9s M %.3f  rozdíl data %s | poměry %s | USp %s | χ² %.1f p %.2e (USp χ² %.0f)' % ('%s_%s' % (F, nm), np.average(Ms, weights=wts), np.round(dat, 3), np.round(pred, 3), np.round(usp, 3), c2, p, out['skupiny']['%s_%s' % (F, nm)]['chi2_usp']), flush=True)
out['verdikt'] = 'O' if all(p > 1e-3 for p in ps) else ('Z' if sum(p < 1e-4 for p in ps) >= 2 else 'N'); out['pmin'] = min(ps)
print('verdikt', out['verdikt'], round(time.time() - t0), 's')
json.dump(out, open('/home/claude/br79/data/fazzari.json', 'w'), indent=1)
