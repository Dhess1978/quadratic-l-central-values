"""34.84 STŘED V · VII — vyhodnocení podle PRAVIDLA_3484.md (zapsána 12. 11. 2026 kolem 19:05 UTC před výpočtem nul D100c).
odv-lu8-rep: vážená hustota samotné D100c (|D| = 10⁸ + 10⁶ … 10⁸ + 2·10⁶, nová D) proti domněnce poměrů (Fazzari 4.8, kusy 34.79);
ODVOZENO, když amplituda sudých a/chyba > 3, (0,80 − a)/chyba > 4 a liché p > 0,001; ZAMÍTNUTO, když a/chyba < 2 nebo (0,80 − a)/chyba < 2.
Popisně: D100 ∪ D100b ∪ D100c a fit bez B10. Postup beze změny proti lu8_analyza.py (34.83): bootstrap 300 a promíchání L 200, chyba = větší z obou,
průmět na tvar odchylky C100 sudých.
Použití: python3 -I lu9_analyza.py --rep SLOZKA_D100c VYSTUP.json
         python3 -I lu9_analyza.py --vse SLOZKA_D100b SLOZKA_D100c VYSTUP.json"""
import sys, os, gzip, json, math, numpy as np, scipy.special as sp
from scipy import stats, optimize
B79, B82 = '/home/claude/br79/data/', '/home/claude/br82/data/'
rng = np.random.default_rng(34841)
REZ = sys.argv[1]; assert REZ in ('--rep', '--vse'); REP = REZ == '--rep'; OUT = sys.argv[-1]
SL = dict(zip(['D100c'] if REP else ['D100b', 'D100c'], sys.argv[2:-1])); JEN = False
EDG = np.linspace(0, 1, 6); BW = 0.2
# ---------------------------------------------------------------- předpověď (stejně jako lu_analyza.py 34.82)
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
# ---------------------------------------------------------------- data: nuly D100 (sada velky) + D100b (Mac 12. 11.), L(½) z PARI (34.72, 34.77)
def cti(P):
    Ds, Gs = [], []
    for fn in sorted(os.listdir(P)):
        if fn.endswith('.jsonl.gz'):
            for line in gzip.open(os.path.join(P, fn), 'rt'): r = json.loads(line); Ds.append(int(r['D'])); Gs.append(r['g'])
    return Ds, Gs
Ds, Gs, zdroj = [], [], []; lv = {}
if not REP:
    Ds, Gs = cti('/home/claude/br63/KS/velky_2026_11_02/D100'); zdroj = ['D100'] * len(Ds)
    f = np.load('/home/claude/br76/data/D100.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
z = np.load('/home/claude/br78/data/D100bc.npz')
for F, P in SL.items():
    Db, Gb = cti(P); Ds += Db; Gs += Gb; zdroj += [F] * len(Db)
    m = z['F'] == F; lv.update(zip(z['D'][m].tolist(), z['L'][m].tolist()))
assert len(set(Ds)) == len(Ds), 'opakovaná D'
bezL = sum(d not in lv for d in Ds); bez_nul = sum(len(g) == 0 for g in Gs)
print('D %d (D100 %d, D100b %d, D100c %d), bez L %d, bez nul %d' % (len(Ds), zdroj.count('D100'), zdroj.count('D100b'), zdroj.count('D100c'), bezL, bez_nul), flush=True)
out = dict(rezim='rep' if REP else 'vse', data=dict(n=len(Ds), D100=zdroj.count('D100'), D100b=zdroj.count('D100b'), D100c=zdroj.count('D100c'), bez_L=bezL, bez_nul=bez_nul, absmin=int(min(abs(d) for d in Ds)), absmax=int(max(abs(d) for d in Ds))), jen_D100=JEN, skup={})
FZ = json.load(open(B79 + 'fazzari.json'))['skupiny']; VN = json.load(open(B79 + 'vaha_nul.json'))
u = np.array(FZ['C100_kl']['rozdil']) - np.array(FZ['C100_kl']['pred'])   # tvar odchylky C100 sudých (34.82)
def amp(r, e): w = u / e ** 2; return float(np.sum(w * r) / np.sum(w * u)), float(1 / math.sqrt(np.sum(u ** 2 / e ** 2)))
for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
    idx = [i for i, d in enumerate(Ds) if (d > 0) == (sg > 0) and d in lv]; n = len(idx)
    L = np.array([lv[Ds[i]] for i in idx]); C = np.zeros((n, 5)); x1 = np.full(n, np.nan)
    for j, i in enumerate(idx):
        t = np.asarray(Gs[i], float)
        if t.size == 0: continue
        x = ((t / 2) * math.log(abs(Ds[i]) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi; C[j] = np.histogram(x, EDG)[0]; x1[j] = x[0]
    st = lambda w: ((w * L) @ C / np.sum(w * L) - w @ C / np.sum(w)) / BW
    d = st(np.ones(n)); eb = np.array([st(rng.poisson(1.0, n).astype(float)) for _ in range(300)]).std(0, ddof=1)
    nul = np.array([((Lp := rng.permutation(L)) @ C / Lp.sum() - C.mean(0)) / BW for _ in range(200)]); en_ = nul.std(0, ddof=1); e = np.maximum(eb, en_)
    l = np.log(np.abs(np.array([Ds[i] for i in idx], float)) / math.pi); pr, nev, M = predikce(l, a)
    c2 = float(np.sum(((d - pr) / e) ** 2)); a_, ea = amp(d - pr, e)
    out['skup'][nm] = dict(n=n, l=float(l.mean()), M=M, L_prumer=float(L.mean()), rozdil=d.tolist(), pred=pr.tolist(), chyba_boot=eb.tolist(), chyba_nul=en_.tolist(), chyba=e.tolist(),
                           chi2=c2, p=float(stats.chi2.sf(c2, 5)), nul_prumer=nul.mean(0).tolist(), nev_dat=(C.mean(0) / BW).tolist(), nev_pred=nev.tolist(), a=a_, e_a=ea, mx1=float(np.nanmean(x1)))
    print('%-3s n %d l %.3f  M %.3f / L̄ %.3f  rozdíl %s  pred %s  χ² %.1f (p %.3g)  amplituda %.3f ± %.3f' % (nm, n, l.mean(), M, L.mean(), np.round(d, 4), np.round(pr, 4), c2, out['skup'][nm]['p'], a_, ea), flush=True)
K8, Z8 = out['skup']['kl'], out['skup']['zap']
zA, zB = K8['a'] / K8['e_a'], (0.80 - K8['a']) / K8['e_a']; out['z_nula'] = zA; out['z_1l'] = zB
if REP: out['verd_rep'] = 'O' if (zA > 3 and zB > 4 and Z8['p'] > 1e-3) else ('Z' if (zA < 2 or zB < 2) else 'N')
print('a/chyba %.2f, (0,80 − a)/chyba %.2f, liché p %.3g%s' % (zA, zB, Z8['p'], ('  → odv-lu8-rep ' + out['verd_rep']) if REP else ''), flush=True)
# ---------------------------------------------------------------- škálování: C100, R10, D100 ∪ D100b, E (bez B10); popisně s B10
E82 = json.load(open(B82 + 'lu.json'))
def bod(F, nm):
    if F == 'D8': g = out['skup'][nm]; return dict(F='D100c' if REP else 'D100+D100b+D100c', l=g['l'], a=g['a'], e=g['e_a'])
    if F == 'E': g = E82['E'][nm]; r = np.array(g['rozdil']) - np.array(g['pred']); e = np.array(g['chyba']); A, eA = amp(r, e); return dict(F='E', l=g['l'], a=A, e=eA)
    g = FZ[F + '_' + nm]; r = np.array(g['rozdil']) - np.array(g['pred']); e = np.maximum(np.array(g['chyba']), np.array(VN[F + '_' + nm]['nul_sd'])); A, eA = amp(r, e); return dict(F=F, l=g['l'], a=A, e=eA)
def fit(B, lC):
    lv_ = np.array([q['l'] for q in B]); av = np.array([q['a'] for q in B]); ev = np.array([q['e'] for q in B])
    fM = lambda x, c, th: c * np.exp(-th * (x - lC)); pM, cM = optimize.curve_fit(fM, lv_, av, p0=[1, 0.3], sigma=ev, absolute_sigma=True, maxfev=20000)
    c2M = float(np.sum(((av - fM(lv_, *pM)) / ev) ** 2)); cL = float(np.sum(av * lC / lv_ / ev ** 2) / np.sum((lC / lv_) ** 2 / ev ** 2)); c2L = float(np.sum(((av - cL * lC / lv_) / ev) ** 2))
    return dict(body=B, M=dict(c=float(pM[0]), theta=float(pM[1]), e_theta=float(math.sqrt(cM[1, 1])), chi2=c2M, p=float(stats.chi2.sf(c2M, len(B) - 2))), L=dict(c=cL, chi2=c2L, p=float(stats.chi2.sf(c2L, len(B) - 1))))
out['skal'] = {}
for nm in ('kl', 'zap'):
    P5 = [bod(F, nm) for F in ('B10', 'C100', 'R10', 'D8', 'E')]; lC = P5[1]['l']
    out['skal'][nm] = dict(bez_B10=fit(P5[1:], lC), vse=fit(P5, lC))
    for k in ('bez_B10', 'vse'):
        S = out['skal'][nm][k]; print('skal %-3s %-7s ' % (nm, k) + ' '.join('%s %.3f±%.3f' % (q['F'], q['a'], q['e']) for q in S['body']) + ' | θ %.3f ± %.3f p %.3g | 1/l χ² %.1f p %.3g' % (S['M']['theta'], S['M']['e_theta'], S['M']['p'], S['L']['chi2'], S['L']['p']), flush=True)
json.dump(out, open(OUT, 'w'), indent=1)
