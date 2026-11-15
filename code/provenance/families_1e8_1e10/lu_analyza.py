"""34.82 STŘED V · V — výpočet tří hypotéz podle PRAVIDLA_3482.md (zapsána před výpočtem).
odv-lu-e10: vážená hustota nízkých nul u |D| ≈ 10¹⁰ (E10 + E100) proti Fazzariho předpovědi (kusy 34.79, nevážená část limitou A → ∞ = vzorec 34.62);
odv-lu-skal: amplituda odchylky sudých B10 … E proti modelům e^{−θ(l − l_C100)} a 1/l;
odv-lu-male: počet L < 0,03 proti a(−3/2)·(2/3)·Res_{−3/2}·ε^{3/2}·n s N z první nuly (N_x1), srovnání s N z prvního momentu.
Vstupy: br82/data/E_nuly.npz, E_L2_*.npz (vlastní AFE), br79/data/{fazzari,fazzari_kusy,fazzari_tridy_kusy,nef}.json/npz, br76 rodiny, br78 h2_vaha (USp tabulka)."""
import json, math, glob, numpy as np, scipy.special as sp, mpmath as mp
from scipy import stats, optimize
B79, B82 = '/home/claude/br79/data/', '/home/claude/br82/data/'
rng = np.random.default_rng(34821)
# ---------------------------------------------------------------- data E: D, L, nuly
EN = np.load(B82 + 'E_nuly.npz'); D = EN['D']; ptr = EN['ptr']; nuly = EN['nuly']; L = np.full(len(D), np.nan)
for f in glob.glob(B82 + 'E_L2_*_2.npz'):
    z = np.load(f); L[z['i']] = z['L']
assert np.all(np.isfinite(L)), 'chybí L'; assert np.all(L > 0), 'L ≤ 0'
EDG = np.linspace(0, 1, 6); BW = 0.2
def xun(t, d, a): return ((t / 2) * np.log(abs(d) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi
K = np.load(B79 + 'fazzari_kusy.npz'); TG = K['TG']; AJ = K['AJ']; NT, NC = len(TG), len(AJ); DT = 0.002
KT = np.load(B79 + 'fazzari_tridy_kusy.npz'); PC = {a: (K['Q%d' % a], K['M%d' % a], KT['U_vse_%d' % a]) for a in (0, 1)}
IMLG = {a: np.imag(sp.loggamma(0.25 + 0.5 * a + 0.5j * TG)) for a in (0, 1)}
def kose(l, a, c):
    tt = np.concatenate([[0], TG]); c0 = 3 * c[0] - 3 * c[1] + c[2]; cc = np.concatenate([[c0], c]); C = np.concatenate([[0], np.cumsum(0.5 * (cc[1:] + cc[:-1]) * DT)])
    th = (0.5 * tt * l + np.concatenate([[0], IMLG[a]])) / math.pi; assert th[-1] > EDG[-1]; te = np.interp(EDG, th, tt); return np.diff(np.interp(te, tt, C))
def predikce(l, a):
    Q, Mk, U = PC[a]; ed = np.linspace(l.min(), l.max() + 1e-9, 21); cnt, _ = np.histogram(l, ed); cen = [l[(l >= ed[i]) & (l < ed[i + 1])].mean() for i in range(20) if cnt[i] > 0]; w_ = cnt[cnt > 0].astype(float)
    num = np.zeros(5); den = 0.0; nu = np.zeros(5); Ms = []
    for w, lc in zip(w_, cen):
        M = np.mean(Mk[:, 0] + np.exp(-AJ * lc) * Mk[:, 1]).real
        E = np.stack([np.ones((NT, NC)), np.exp(-AJ * lc)[None, :].repeat(NT, 0), np.exp(-1j * TG * lc)[:, None].repeat(NC, 1), np.exp(-(AJ[None, :] + 1j * TG[:, None]) * lc)], -1)
        cw = np.real(np.mean(np.sum(Q * E, -1), 1)) / math.pi; cu = np.real(U[:, 0] + np.exp(-1j * TG * lc) * U[:, 1]) / math.pi
        Ms.append(M); num += w * (M * BW + kose(lc, a, cw)); den += w * M * BW; nu += w * (BW + kose(lc, a, cu))
    return num / den - nu / (w_.sum() * BW), nu / (w_.sum() * BW), float(np.average(Ms, weights=w_))
out = {'E': {}}
H2 = json.load(open('/home/claude/br78/data/h2_vaha.json')); U = {int(k): v[0] for k, v in H2['usp'].items()}
def N_z_x1(mx):   # inverze E[x₁](N) USp(2N), lineárně mezi celými N
    Ns = sorted(U); v = [U[n] for n in Ns]; return float(np.interp(mx, v, Ns))
for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
    m = np.nonzero((D > 0) == (sg > 0))[0]; n = len(m); C = np.zeros((n, 5)); x1 = np.zeros(n)
    for j, i in enumerate(m):
        x = xun(nuly[ptr[i]:ptr[i + 1]], D[i], a); C[j] = np.histogram(x, EDG)[0]; x1[j] = x[0]
    Lm = L[m]; st = lambda w: ((w * Lm) @ C / np.sum(w * Lm) - w @ C / np.sum(w)) / BW
    d = st(np.ones(n)); Bt = np.array([st(rng.poisson(1.0, n).astype(float)) for _ in range(300)]); eb = Bt.std(0, ddof=1)
    nul = np.array([((Lp := rng.permutation(Lm)) @ C / Lp.sum() - C.mean(0)) / BW for _ in range(200)]); en_ = nul.std(0, ddof=1); e = np.maximum(eb, en_)
    l = np.log(np.abs(D[m]) / math.pi); pr, nev, M = predikce(l, a); c2 = float(np.sum(((d - pr) / e) ** 2)); p = float(stats.chi2.sf(c2, 5))
    c2b = float(np.sum(((d - pr) / eb) ** 2))
    out['E'][nm] = dict(n=n, l=float(l.mean()), M=M, L_prumer=float(Lm.mean()), rozdil=d.tolist(), pred=pr.tolist(), chyba_boot=eb.tolist(), chyba_nul=en_.tolist(), chyba=e.tolist(), chi2=c2, p=p, chi2_boot=c2b, p_boot=float(stats.chi2.sf(c2b, 5)),
                        nev_dat=(C.mean(0) / BW).tolist(), nev_pred=nev.tolist(), nul_prumer=nul.mean(0).tolist(), mx1=float(x1.mean()), se_x1=float(x1.std(ddof=1) / math.sqrt(n)), N_x1=N_z_x1(float(x1.mean())),
                        pod003=int((Lm < 0.03).sum()), pod001=int((Lm < 0.01).sum()), pod01=int((Lm < 0.1).sum()), Lmin=float(Lm.min()))
    print('E %-3s n %d l %.2f M %.3f / L̄ %.3f  rozdíl %s  pred %s  χ² %.1f (p %.2e; boot %.1f)  N_x1 %.2f' % (nm, n, l.mean(), M, Lm.mean(), np.round(d, 3), np.round(pr, 3), c2, p, c2b, out['E'][nm]['N_x1']), flush=True)
pv = [out['E'][s]['p'] for s in ('kl', 'zap')]; out['verd_e10'] = 'O' if all(p > 1e-3 for p in pv) else ('Z' if any(p < 1e-4 for p in pv) else 'N')
# ---------------------------------------------------------------- odv-lu-skal: amplitudy
FZ = json.load(open(B79 + 'fazzari.json'))['skupiny']; VN = json.load(open(B79 + 'vaha_nul.json')); out['skal'] = {}
for nm in ('kl', 'zap'):
    u = np.array(FZ['C100_kl']['rozdil']) - np.array(FZ['C100_kl']['pred'])   # tvar odchylky C100 sudých (pro obě znaménka týž)
    pts = []
    for F in ('B10', 'C100', 'R10', 'D100', 'E'):
        if F == 'E': g = out['E'][nm]; r = np.array(g['rozdil']) - np.array(g['pred']); e = np.array(g['chyba']); l = g['l']
        else: g = FZ[F + '_' + nm]; r = np.array(g['rozdil']) - np.array(g['pred']); e = np.maximum(np.array(g['chyba']), np.array(VN[F + '_' + nm]['nul_sd'])); l = g['l']
        w = u / e ** 2; a_ = float(np.sum(w * r) / np.sum(w * u)); ea = float(1 / math.sqrt(np.sum(u ** 2 / e ** 2))); pts.append(dict(F=F, l=l, a=a_, e=ea))
    lv = np.array([q['l'] for q in pts]); av = np.array([q['a'] for q in pts]); ev = np.array([q['e'] for q in pts]); lC = pts[1]['l']
    fM = lambda l, c, th: c * np.exp(-th * (l - lC)); pM, cM = optimize.curve_fit(fM, lv, av, p0=[1, 0.3], sigma=ev, absolute_sigma=True, maxfev=20000)
    c2M = float(np.sum(((av - fM(lv, *pM)) / ev) ** 2)); cL = float(np.sum(av * lC / lv / ev ** 2) / np.sum((lC / lv) ** 2 / ev ** 2)); c2L = float(np.sum(((av - cL * lC / lv) / ev) ** 2))
    out['skal'][nm] = dict(body=pts, M=dict(c=float(pM[0]), theta=float(pM[1]), e_theta=float(math.sqrt(cM[1, 1])), chi2=c2M, p=float(stats.chi2.sf(c2M, 3))), L=dict(c=cL, chi2=c2L, p=float(stats.chi2.sf(c2L, 4))))
    print('skal %-3s ' % nm + ' '.join('%s %.3f±%.3f' % (q['F'], q['a'], q['e']) for q in pts) + ' | M θ %.3f ± %.3f χ² %.1f p %.3g | L χ² %.1f p %.3g' % (pM[1], math.sqrt(cM[1, 1]), c2M, out['skal'][nm]['M']['p'], c2L, out['skal'][nm]['L']['p']), flush=True)
S = out['skal']['kl']; out['verd_skal'] = 'O' if (S['M']['p'] > 0.01 and S['L']['p'] < 1e-3) else ('Z' if S['L']['p'] > 0.01 else 'N')
# ---------------------------------------------------------------- odv-lu-male
LH = json.load(open('/home/claude/br76/data/lh_ref.json')); a32 = LH['male']['a32'] if 'male' in LH else None
mp.mp.dps = 25; lnG = lambda z: float(mp.log(mp.barnesg(z)))
lnres = lambda N: -3 * N * math.log(2) + lnG(2 * N + 2) - lnG(N + 2) - lnG(N + 1.5) + lnG(1.5) + lnG(N) - lnG(2 * N + 0.5) + lnG(N + 0.5)
K1 = LH['k1']; a1 = LH['a1']
lnE = lambda N: 2 * N * math.log(2) + math.lgamma(N + 1.5) - math.lgamma(1.5) - math.lgamma(2 * N + 2) + math.lgamma(N + 2)
def Nk1(mq):
    t = math.log(mq / a1); lo, hi = 0.01, 60.0
    for _ in range(200):
        mm = (lo + hi) / 2; lo, hi = (mm, hi) if lnE(mm) < t else (lo, mm)
    return (lo + hi) / 2
if a32 is None:
    s = open('/home/claude/br76/lh3476data.js').read(); a32 = json.loads(s[s.index('{'):s.rindex('}') + 1])['male']['a32']
NF = json.load(open(B79 + 'nef.json'))['skupiny']; EPS = 0.03; out['male'] = {'a32': a32, 'eps': EPS, 'skupiny': {}}; zs = []
for F in ('B10', 'C100', 'R10', 'D100', 'E'):
    if F == 'E': ff = {'D': D, 'L': L}
    else: z = np.load('/home/claude/br76/data/' + F + '.npz'); ff = {'D': z['D'].astype(np.int64), 'L': z['L']}
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        m = (ff['D'] > 0) == (sg > 0); Lm = ff['L'][m]; n = int(m.sum()); obs = {e: int((Lm < e).sum()) for e in (0.1, 0.03, 0.01)}
        Nx = out['E'][nm]['N_x1'] if F == 'E' else NF[F + '_' + nm]['N_x1']
        l = np.log(np.abs(ff['D'][m]).astype(float)); Nk = Nk1(K1['c0'][a] + K1['c1'] * float(l.mean()))
        Cx = a32 * math.exp(lnres(Nx)) * 2 / 3; Ck = a32 * math.exp(lnres(Nk)) * 2 / 3
        px = {e: Cx * e ** 1.5 * n for e in (0.1, 0.03, 0.01)}; pk = {e: Ck * e ** 1.5 * n for e in (0.1, 0.03, 0.01)}
        z_ = (obs[EPS] - px[EPS]) / math.sqrt(px[EPS]); zs.append(z_)
        out['male']['skupiny'][F + '_' + nm] = dict(n=n, N_x1=Nx, N_k1=Nk, C_x1=Cx, C_k1=Ck, obs={str(e): v for e, v in obs.items()}, pred_x1={str(e): v for e, v in px.items()}, pred_k1={str(e): v for e, v in pk.items()}, z=z_, z_k1=(obs[EPS] - pk[EPS]) / math.sqrt(pk[EPS]))
        print('male %-9s n %7d N_x1 %.2f N_k1 %.2f  pod 0,03: %d  pred(N_x1) %.1f  pred(N_k1) %.1f  z %+.1f (z_k1 %+.1f)' % (F + '_' + nm, n, Nx, Nk, obs[EPS], px[EPS], pk[EPS], z_, out['male']['skupiny'][F + '_' + nm]['z_k1']), flush=True)
zm = max(abs(z) for z in zs); out['male']['zmax'] = zm; out['verd_male'] = 'O' if zm < 3 else ('Z' if zm >= 4 else 'N')
print('verdikty e10 %s skal %s male %s' % (out['verd_e10'], out['verd_skal'], out['verd_male']))
json.dump(out, open(B82 + 'lu.json', 'w'), indent=1)
