"""34.79 popisně k odv-lr-vaha (David 12. 11. se ptal, zda odchylka vážené hustoty není chybou přerušeného kontrolního balíku). Předpověď domněnky poměrů s vahou L(½) zvlášť pro třídy D podle χ_D(2):
D ≡ 1 mod 8 (χ(2) = +1), D ≡ 5 mod 8 (χ(2) = −1), 4 | D (χ(2) = 0) — mění se jen lokální faktor u p = 2; třída 'vse' = celá rodina (kontrola proti fazzari.json).
Nevážená hustota týmž kódem limitou A → ∞ (L(½ + A) → 1) — nezávislá kontrola vzorce 34.62. Data: rozdíl vážené − nevážené po třídách přímo z balíků nul, bootstrap 300×.
Hodnoty ζ se počítají jednou a ukládají (obnovitelné po restartu)."""
import os, sys, gzip, json, math, csv, time, numpy as np, mpmath as mp, scipy.special as sp
from scipy import stats
t0 = time.time(); mp.mp.dps = 30; OUT = '/home/claude/br79/data/'
Mp = 200000; s = np.ones(Mp + 1, bool); s[:2] = False
for p in range(2, int(Mp ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0].astype(float); lp = np.log(P); RP = P / (P + 1); assert P[0] == 2
pw = lambda e: np.exp(-e * lp)
def br_fam(x, y, w, eps):
    """lokální faktor (závislý na χ(p)) a jeho derivace podle C; eps None = průměr rodiny, jinak pevná hodnota χ(2)"""
    pp, mm = (1 - x) * (1 - y), (1 + x) * (1 + y); E0 = 0.5 * (1 / pp + 1 / mm) - 1; O = 0.5 * (1 / pp - 1 / mm); br = 1 + RP * (E0 - w * O)
    dE0 = 0.5 * (1 / ((1 - x) * (1 - y) ** 2) - 1 / ((1 + x) * (1 + y) ** 2)); dO = 0.5 * (1 / ((1 - x) * (1 - y) ** 2) + 1 / ((1 + x) * (1 + y) ** 2)); dbr = RP * (dE0 - w * dO) * (-lp * y)
    if eps is not None:
        br = br.copy(); dbr = dbr.copy(); X, Y, W, L2 = x[0], y[0], w[0], lp[0]
        if eps == 1: br[0] = (1 - W) / ((1 - X) * (1 - Y)); dbr[0] = (1 - W) / ((1 - X) * (1 - Y) ** 2) * (-L2 * Y)
        elif eps == -1: br[0] = (1 + W) / ((1 + X) * (1 + Y)); dbr[0] = (1 + W) / ((1 + X) * (1 + Y) ** 2) * (L2 * Y)
        else: br[0] = 1.0; dbr[0] = 0.0
    return br, dbr
def A2(A, C, D, eps, inf=False):
    y, w = pw(0.5 + C), pw(0.5 + D); c2, cd = pw(1 + 2 * C), pw(1 + C + D)
    if inf: x = np.zeros_like(y); f1 = np.log(1 - c2) - np.log(1 - cd); d1 = 2 * lp * c2 / (1 - c2) - lp * cd / (1 - cd)
    else:
        x = pw(0.5 + A); a2, ac, ad = pw(1 + 2 * A), pw(1 + A + C), pw(1 + A + D)
        f1 = np.log(1 - a2) + np.log(1 - c2) + np.log(1 - ac) - np.log(1 - ad) - np.log(1 - cd); d1 = 2 * lp * c2 / (1 - c2) + lp * ac / (1 - ac) - lp * cd / (1 - cd)
    br, dbr = br_fam(x, y, w, eps); val = np.exp(np.sum(f1 + np.log(br))); return val, val * np.sum(d1 + dbr / br)
def A1(A, eps):
    x = pw(0.5 + A); E0 = 0.5 * (1 / (1 - x) + 1 / (1 + x)) - 1; f = 1 + RP * E0
    if eps is not None: f = f.copy(); f[0] = {1: 1 / (1 - x[0]), -1: 1 / (1 + x[0]), 0: 1.0}[eps]
    return np.exp(np.sum(np.log(1 - pw(1 + 2 * A)) + np.log(f)))
Z = lambda z: complex(mp.zeta(mp.mpc(z.real, z.imag))); Zd = lambda z: complex(mp.zeta(mp.mpc(z.real, z.imag), derivative=1))
def gS(z, a): return complex(mp.gamma(mp.mpc((1 - z.real + a) / 2, -z.imag / 2)) / mp.gamma(mp.mpc((z.real + a) / 2, z.imag / 2)))
NC, H = 16, 0.0037; AJ = np.array([H * complex(math.cos(2 * math.pi * (j + 0.5) / NC), math.sin(2 * math.pi * (j + 0.5) / NC)) for j in range(NC)])
DT = 0.002; TG = np.round(np.arange(DT, 1.40 + 1e-9, DT), 4); NT = len(TG)
# ---- cache ζ a Γ (obnovitelné)
CF = OUT + 'fazzari_tridy_zeta.npz'
if os.path.exists(CF): Cz = dict(np.load(CF)); print('cache ζ načtena', flush=True)
else:
    Cz = {k: np.zeros((NT, NC), complex) for k in ('zpd', 'zdpd', 'zmd', 'zdmd', 'zpmd', 'zmmd', 'g0', 'g1')}
    Cz.update({k: np.zeros(NT, complex) for k in ('zd2', 'zdd2', 'zm2', 'gd0', 'gd1')}); Cz.update({k: np.zeros(NC, complex) for k in ('za', 'zma', 'ga0', 'ga1')})
    for j, A in enumerate(AJ): Cz['za'][j] = Z(1 + 2 * A); Cz['zma'][j] = Z(1 - 2 * A); Cz['ga0'][j] = gS(0.5 + A, 0); Cz['ga1'][j] = gS(0.5 + A, 1)
    for i, t in enumerate(TG):
        D = 1j * float(t); Cz['zd2'][i] = Z(1 + 2 * D); Cz['zdd2'][i] = Zd(1 + 2 * D); Cz['zm2'][i] = Z(1 - 2 * D); Cz['gd0'][i] = gS(0.5 + D, 0); Cz['gd1'][i] = gS(0.5 + D, 1)
        for j, A in enumerate(AJ):
            Cz['zpd'][i, j] = Z(1 + A + D); Cz['zdpd'][i, j] = Zd(1 + A + D); Cz['zmd'][i, j] = Z(1 - A + D); Cz['zdmd'][i, j] = Zd(1 - A + D); Cz['zpmd'][i, j] = Z(1 + A - D); Cz['zmmd'][i, j] = Z(1 - A - D)
            Cz['g0'][i, j] = gS(0.5 + A + D, 0); Cz['g1'][i, j] = gS(0.5 + A + D, 1)
        if i % 100 == 0: print('ζ t', t, round(time.time() - t0), 's', flush=True)
    np.savez(CF, **Cz)
# ---- kusy pro třídy: A2 nezávisí na a ani na třídě kromě faktoru u p = 2 → jeden průchod přes prvočísla, třídy a parity skalární opravou u p = 2
TR = {'vse': None, '1mod8': 1, '5mod8': -1, '4deli': 0}
def A2zakl(A, C, D, inf=False):
    """val a Σ(d1 + dbr/br) pro rodinu + veličiny u p = 2 (x, y, w, br, dbr)"""
    y, w = pw(0.5 + C), pw(0.5 + D); c2, cd = pw(1 + 2 * C), pw(1 + C + D)
    if inf: x = np.zeros_like(y); f1 = np.log(1 - c2) - np.log(1 - cd); d1 = 2 * lp * c2 / (1 - c2) - lp * cd / (1 - cd)
    else:
        x = pw(0.5 + A); a2, ac, ad = pw(1 + 2 * A), pw(1 + A + C), pw(1 + A + D)
        f1 = np.log(1 - a2) + np.log(1 - c2) + np.log(1 - ac) - np.log(1 - ad) - np.log(1 - cd); d1 = 2 * lp * c2 / (1 - c2) + lp * ac / (1 - ac) - lp * cd / (1 - cd)
    br, dbr = br_fam(x, y, w, None); val = np.exp(np.sum(f1 + np.log(br))); return np.array([val, np.sum(d1 + dbr / br), x[0], y[0], w[0], br[0], dbr[0]])
def A2tr(z, eps):
    val, S, X, Y, W, bf, dbf = z
    if eps is None: return val, val * S
    if eps == 1: bc = (1 - W) / ((1 - X) * (1 - Y)); dbc = (1 - W) / ((1 - X) * (1 - Y) ** 2) * (-lp[0] * Y)
    elif eps == -1: bc = (1 + W) / ((1 + X) * (1 + Y)); dbc = (1 + W) / ((1 + X) * (1 + Y) ** 2) * (lp[0] * Y)
    else: bc, dbc = 1.0, 0.0
    v = val * bc / bf; return v, v * (S - dbf / bf + dbc / bc)
ZF = OUT + 'fazzari_tridy_A2.npz'
if os.path.exists(ZF): Z2 = np.load(ZF)['Z2']; ZI = np.load(ZF)['ZI']; print('cache A2 načtena', flush=True)
else:
    Z2 = np.zeros((NT, NC, 4, 7), complex); ZI = np.zeros((NT, 2, 7), complex)
    for i, t in enumerate(TG):
        D = 1j * float(t); ZI[i, 0] = A2zakl(0, D, D, True); ZI[i, 1] = A2zakl(0, -D, D, True)
        for j, A in enumerate(AJ): Z2[i, j] = [A2zakl(A, D, D), A2zakl(-A, D, D), A2zakl(A, -D, D), A2zakl(-A, -D, D)]
        if i % 100 == 0: print('A2 t', t, round(time.time() - t0), 's', flush=True)
    np.savez(ZF, Z2=Z2, ZI=ZI)
def kusy(eps, a):
    ga, gd, g2 = Cz['ga%d' % a], Cz['gd%d' % a], Cz['g%d' % a]; za, zma = Cz['za'], Cz['zma']
    Mk = np.array([[za[j] * A1(A, eps), ga[j] * zma[j] * A1(-A, eps)] for j, A in enumerate(AJ)]); Q = np.zeros((NT, NC, 4), complex); U = np.zeros((NT, 2), complex)
    for i, t in enumerate(TG):
        lz = Cz['zdd2'][i] / Cz['zd2'][i]; zm2 = Cz['zm2'][i]
        v, dv = A2tr(ZI[i, 0], eps); v3, _ = A2tr(ZI[i, 1], eps); U[i] = [v * lz + dv, gd[i] * (-v3 * zm2)]
        for j in range(NC):
            vA, dA = A2tr(Z2[i, j, 0], eps); vmA, dmA = A2tr(Z2[i, j, 1], eps); v3, _ = A2tr(Z2[i, j, 2], eps); v4, _ = A2tr(Z2[i, j, 3], eps)
            zpd, zmd = Cz['zpd'][i, j], Cz['zmd'][i, j]
            q1 = vA * za[j] / zpd * (lz * zpd + Cz['zdpd'][i, j]) + dA * za[j]; q2 = vmA * zma[j] / zmd * (lz * zmd + Cz['zdmd'][i, j]) + dmA * zma[j]
            q3 = -v3 * za[j] * zm2 * Cz['zpmd'][i, j] / zpd; q4 = -v4 * zma[j] * zm2 * Cz['zmmd'][i, j] / zmd
            Q[i, j] = [q1, ga[j] * q2, gd[i] * q3, g2[i, j] * q4]
    return Q, Mk, U
# kontrola opravy u p = 2 proti přímému výpočtu (3 body)
for eps in (1, -1, 0):
    A, D = AJ[3], 0.37j; v1 = A2(A, D, D, eps); v2 = A2tr(A2zakl(A, D, D), eps); assert abs(v1[0] - v2[0]) < 1e-12 * abs(v1[0]) and abs(v1[1] - v2[1]) < 1e-10 * abs(v1[1]), (eps, v1, v2)
    v1 = A2(0, -D, D, eps, inf=True); v2 = A2tr(A2zakl(0, -D, D, True), eps); assert abs(v1[0] - v2[0]) < 1e-12 * abs(v1[0]), eps
print('oprava u p = 2 = přímý výpočet', flush=True)
IMLG = {k: np.imag(sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
def kose(l, a, c, edges):
    tt = np.concatenate([[0], TG]); c0 = 3 * c[0] - 3 * c[1] + c[2]; cc = np.concatenate([[c0], c]); C = np.concatenate([[0], np.cumsum(0.5 * (cc[1:] + cc[:-1]) * DT)])
    th = (0.5 * tt * l + np.concatenate([[0], IMLG[a]])) / math.pi; te = np.interp(edges, th, tt); return np.diff(np.interp(te, tt, C))
EDG = np.linspace(0, 1, 6); BW = 0.2
def predikce(K, a, l):
    Q, Mk, U = K; ed = np.linspace(l.min(), l.max() + 1e-9, 21); cnt, _ = np.histogram(l, ed); cen = [l[(l >= ed[i]) & (l < ed[i + 1])].mean() for i in range(20) if cnt[i] > 0]; wts = cnt[cnt > 0].astype(float)
    num = np.zeros(5); den = 0.0; nu = np.zeros(5); Ms = []
    for w, lc in zip(wts, cen):
        M = np.mean(Mk[:, 0] + np.exp(-AJ * lc) * Mk[:, 1]).real
        E = np.stack([np.ones((NT, NC)), np.exp(-AJ * lc)[None, :].repeat(NT, 0), np.exp(-1j * TG * lc)[:, None].repeat(NC, 1), np.exp(-(AJ[None, :] + 1j * TG[:, None]) * lc)], -1)
        cw = np.real(np.mean(np.sum(Q * E, -1), 1)) / math.pi; cu = np.real(U[:, 0] + np.exp(-1j * TG * lc) * U[:, 1]) / math.pi
        Ms.append(M); num += w * (M * BW + kose(lc, a, cw, EDG)); den += w * M * BW; nu += w * (BW + kose(lc, a, cu, EDG))
    vaz = num / den; nev = nu / (wts.sum() * BW); return vaz - nev, vaz, nev, float(np.average(Ms, weights=wts))
KF = OUT + 'fazzari_tridy_kusy.npz'
if os.path.exists(KF): KK = dict(np.load(KF)); KS = {k: (KK['Q_' + k], KK['M_' + k], KK['U_' + k]) for k in [x[2:] for x in KK if x.startswith('Q_')]}
else: KS = {}
for tn, eps in TR.items():
    for a in (0, 1):
        k = '%s_%d' % (tn, a)
        if k in KS: continue
        KS[k] = kusy(eps, a); print('kusy', k, round(time.time() - t0), 's', flush=True)
        np.savez(KF, **{p + kk: v[i] for kk, v in KS.items() for i, p in enumerate(('Q_', 'M_', 'U_'))})
# ---- data po třídách
SRC = {'B10': '/home/claude/br63/KS/laborator_male_2026_11_02/B10', 'C100': '/home/claude/br63/KS/velky_2026_11_02/C100', 'R10': '/home/claude/br63/KS/velky_2026_11_02/R10', 'D100': '/home/claude/br63/KS/velky_2026_11_02/D100'}
FZ = json.load(open(OUT + 'fazzari.json'))['skupiny']; rng = np.random.default_rng(347912); out = {'tridy': list(TR), 'skupiny': {}}
for F, Pth in SRC.items():
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist())); Ds, Gs = [], []
    for fn in sorted(os.listdir(Pth)):
        if fn.endswith('.jsonl.gz'):
            for line in gzip.open(os.path.join(Pth, fn), 'rt'): r = json.loads(line); Ds.append(int(r['D'])); Gs.append(r['g'])
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        idx = [i for i, d in enumerate(Ds) if (d > 0) == (sg > 0) and d in lv]; D = np.array([Ds[i] for i in idx], np.int64); L = np.array([lv[Ds[i]] for i in idx]); C = np.zeros((len(idx), 5))
        for j, i in enumerate(idx):
            t = np.asarray(Gs[i], float); x = ((t / 2) * math.log(abs(Ds[i]) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi; C[j] = np.histogram(x, EDG)[0]
        key = F + '_' + nm; out['skupiny'][key] = {}
        for tn in TR:
            m = np.ones(len(D), bool) if tn == 'vse' else ((D % 8 == 1) if tn == '1mod8' else (D % 8 == 5) if tn == '5mod8' else (D % 4 == 0))
            Cm, Lm = C[m], L[m]; st = lambda w: ((w * Lm) @ Cm / np.sum(w * Lm) - w @ Cm / np.sum(w)) / BW
            d = st(np.ones(m.sum())); B = np.array([st(rng.poisson(1.0, m.sum()).astype(float)) for _ in range(300)]); e = B.std(0, ddof=1)
            nev_dat = Cm.sum(0) / m.sum() / BW; l = np.log(np.abs(D[m]) / math.pi); pr, vz, nv, M = predikce(KS['%s_%d' % (tn, a)], a, l)
            c2 = float(np.sum(((d - pr) / e) ** 2)); out['skupiny'][key][tn] = dict(n=int(m.sum()), M=M, L_prumer=float(Lm.mean()), rozdil=d.tolist(), chyba=e.tolist(), pred=pr.tolist(), vaz_pred=vz.tolist(), nev_pred=nv.tolist(), nev_dat=nev_dat.tolist(), chi2=c2, p=float(stats.chi2.sf(c2, 5)))
        o = out['skupiny'][key]; print('%-9s ' % key + ' | '.join('%s n %d M %.2f/%.2f χ² %.1f' % (tn, o[tn]['n'], o[tn]['M'], o[tn]['L_prumer'], o[tn]['chi2']) for tn in TR) + ' | vse−fazzari %.1e, nev−34.62 %.1e' % (np.max(np.abs(np.array(o['vse']['pred']) - np.array(FZ[key]['pred']))), np.max(np.abs(np.array(o['vse']['nev_pred']) - np.array(FZ[key]['nevazena_pred'])))), flush=True)
json.dump(out, open(OUT + 'fazzari_tridy.json', 'w'), indent=1); print('hotovo', round(time.time() - t0), 's')
