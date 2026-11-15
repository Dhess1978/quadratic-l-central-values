"""34.85a STŘED V · VIII — rodina S7 (7,5·10⁶ < |D| ≤ 10⁷, nuly z Macu) podle PRAVIDLA_3485.md: část B (zapsaná před výpočtem) a doplněk D
(zapsaný po výsledku odv-lv-nez a PŘED daty S7).
odv-lv-s7k: opravená předpověď domněnky poměrů s vahou L(½) popisuje S7 pro obě znaménka — ODVOZENO, když p(χ², 5 košů) > 0,001 u sudých i lichých;
ZAMÍTNUTO, když p < 10⁻⁴ u kteréhokoli znaménka; jinak NELZE ODVODIT.
odv-lv-s7 a odv-lv-tvar (část B) se vyhodnotí tak, jak jsou zapsané, ale s opravenou předpovědí (čekaný výsledek ZAMÍTNUTO, předpoklad — člen sudých — padl).
Postup jako 34.83/34.84 (lu9_analyza.py): rozdíl vážené a nevážené hustoty v 5 koších x ∈ [0, 1]; chyba = větší z (bootstrap 300×, promíchání L 200×);
amplituda = vážený průmět na tvar dřívější odchylky C100 sudých. Předpověď: nezávislá implementace (nez_porovnani.predikce, 40 košů l).
Použití: python3 -I s7_analyza.py RODINA SLOZKA_NUL HODNOTY.npz[,…] VYSTUP.json"""
import sys, os, gzip, json, math, time, numpy as np, scipy.special as sp
from scipy import stats
sys.path.insert(0, '/home/claude/br85/py'); import nez_porovnani as NP
t0 = time.time(); F, ZD, HV, OUT = sys.argv[1:5]; rng = np.random.default_rng(34852)
EDG = np.linspace(0, 1, 6); BW = 0.2
# ---------------------------------------------------------------- chybná předpověď 34.79–34.84 (pro záznam a pro obraz „chybná × opravená“), beze změny podle lu9_analyza.py
B79 = '/home/claude/br79/data/'; K0 = np.load(B79 + 'fazzari_kusy.npz'); TG0 = K0['TG']; AJ0 = K0['AJ']; NT0, NC0 = len(TG0), len(AJ0); DT0 = 0.002
KT0 = np.load(B79 + 'fazzari_tridy_kusy.npz'); PC0 = {a: (K0['Q%d' % a], K0['M%d' % a], KT0['U_vse_%d' % a]) for a in (0, 1)}
IMLG0 = {a: np.imag(sp.loggamma(0.25 + 0.5 * a + 0.5j * TG0)) for a in (0, 1)}
def kose0(l, a, c):
    tt = np.concatenate([[0], TG0]); c0 = 3 * c[0] - 3 * c[1] + c[2]; cc = np.concatenate([[c0], c]); C = np.concatenate([[0], np.cumsum(0.5 * (cc[1:] + cc[:-1]) * DT0)])
    th = (0.5 * tt * l + np.concatenate([[0], IMLG0[a]])) / math.pi; assert th[-1] > EDG[-1]; te = np.interp(EDG, th, tt); return np.diff(np.interp(te, tt, C))
def predikce_chybna(l, a):
    Q, Mk, U = PC0[a]; ed = np.linspace(l.min(), l.max() + 1e-9, 21); cnt, _ = np.histogram(l, ed); cen = [l[(l >= ed[i]) & (l < ed[i + 1])].mean() for i in range(20) if cnt[i] > 0]; w_ = cnt[cnt > 0].astype(float)
    num = np.zeros(5); den = 0.0; nu = np.zeros(5)
    for w, lc in zip(w_, cen):
        M = np.mean(Mk[:, 0] + np.exp(-AJ0 * lc) * Mk[:, 1]).real
        E = np.stack([np.ones((NT0, NC0)), np.exp(-AJ0 * lc)[None, :].repeat(NT0, 0), np.exp(-1j * TG0 * lc)[:, None].repeat(NC0, 1), np.exp(-(AJ0[None, :] + 1j * TG0[:, None]) * lc)], -1)
        cw = np.real(np.mean(np.sum(Q * E, -1), 1)) / math.pi; cu = np.real(U[:, 0] + np.exp(-1j * TG0 * lc) * U[:, 1]) / math.pi
        num += w * (M * BW + kose0(lc, a, cw)); den += w * M * BW; nu += w * (BW + kose0(lc, a, cu))
    return num / den - nu / (w_.sum() * BW)
# ---------------------------------------------------------------- nuly
Dl, cnt, gl = [], [], []
for fn in sorted(os.listdir(ZD)):
    if not fn.endswith('.jsonl.gz'): continue
    for line in gzip.open(os.path.join(ZD, fn), 'rt'):
        r = json.loads(line); Dl.append(int(r['D'])); g = r['g']; cnt.append(len(g)); gl.extend(g)
D = np.array(Dl, np.int64); cnt = np.array(cnt); G = np.array(gl, float); own = np.repeat(np.arange(len(D)), cnt)
assert len(np.unique(D)) == len(D), 'opakovaná D'
print(F, 'D', len(D), 'nul', len(G), 'bez nul', int((cnt == 0).sum()), round(time.time() - t0), 's', flush=True)
# ---------------------------------------------------------------- hodnoty L(½) (PARI)
HD, HL = [], []
for f in HV.split(','):
    z = np.load(f); HD.append(z['D']); HL.append(z['L'])
HD = np.concatenate(HD); HL = np.concatenate(HL); s = np.argsort(HD); HD, HL = HD[s], HL[s]
ix = np.clip(np.searchsorted(HD, D), 0, len(HD) - 1); maL = HD[ix] == D; L = np.where(maL, HL[ix], np.nan)
out = dict(rodina=F, data=dict(n=int(len(D)), nul=int(len(G)), bez_L=int((~maL).sum()), bez_nul=int((cnt == 0).sum()), absmin=int(np.abs(D).min()), absmax=int(np.abs(D).max()),
                               L_min=float(np.nanmin(L)), L_nekladne=int((L[maL] <= 0).sum())), skup={})
if F == 'S7':      # úplnost: všechny základní D s 7,5·10⁶ < |D| ≤ 10⁷ z tabulky hodnot
    cil = HD[(np.abs(HD) > 7500000) & (np.abs(HD) <= 10000000)]; out['data']['cil'] = int(len(cil)); out['data']['chybi'] = int(len(np.setdiff1d(cil, D))); out['data']['navic'] = int(len(np.setdiff1d(D, cil)))
    print('úplnost: cíl %d, chybí %d, navíc %d' % (len(cil), out['data']['chybi'], out['data']['navic']), flush=True)
    # počty různých D s nulami (rodina R10 je náhodný výběr z 10⁶ < |D| ≤ 10⁷ a překrývá se s C100 a S7)
    CV = '/home/claude/pub5/repo/data/central_values/'; Vv = np.concatenate([np.load(CV + f)['D'] for f in sorted(os.listdir(CV)) if f.startswith('V_part')])
    Bm = Vv[(np.abs(Vv) > 10 ** 4) & (np.abs(Vv) <= 10 ** 5)]; Cm = Vv[(np.abs(Vv) > 10 ** 6) & (np.abs(Vv) <= 2 * 10 ** 6)]; Rm = np.load('/home/claude/br76/data/R10.npz')['D']
    Om = [np.load(CV + n + '.npz')['D'] for n in ('D100', 'D100b', 'D100c', 'E10_E100')]; U0 = np.unique(np.concatenate([Bm, Cm, Rm] + Om))
    out['prekryv'] = dict(R10_C100=int(len(np.intersect1d(Rm, Cm))), R10_S7=int(len(np.intersect1d(Rm, D))), ruznych_bez_S7=int(len(U0)), ruznych=int(len(np.unique(np.concatenate([U0, D])))),
                          po_rodinach_bez_S7=int(len(Bm) + len(Cm) + len(Rm) + sum(len(x) for x in Om)))
    print('překryv', out['prekryv'], flush=True)
R = json.load(open('/home/claude/pub5/repo/data/results/study_v_results.json'))['weighted_density']
u = np.array(R['amplitude_shape_u'])      # tvar dřívější odchylky C100 sudých (data − chybná předpověď)
eu = np.array(R['corrected']['groups']['C100_even']['error'])      # chyba dat C100 sudých (pro odv-lv-tvar „chyby obou“)
def amp(r, e): w = u / e ** 2; return float(np.sum(w * r) / np.sum(w * u)), float(1 / math.sqrt(np.sum(u ** 2 / e ** 2)))
for nm, sg, a in (('even', 1, 0), ('odd', -1, 1)):
    m = ((D > 0) == (sg > 0)) & maL; n = int(m.sum()); pos = np.full(len(D), -1); pos[m] = np.arange(n)
    Dm = D[m].astype(float); Lm = L[m]; l = np.log(np.abs(Dm) / math.pi)
    k = pos[own]; ok = k >= 0; t = G[ok]; k = k[ok]
    x = ((t / 2) * l[k] + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi
    b = np.floor(x / BW).astype(int); inb = (x >= 0) & (b < 5); C = np.zeros((n, 5)); np.add.at(C, (k[inb], b[inb]), 1.0)
    first = np.r_[True, k[1:] != k[:-1]]; x1 = float(np.mean(x[first]))
    st = lambda w: ((w * Lm) @ C / np.sum(w * Lm) - w @ C / np.sum(w)) / BW
    d = st(np.ones(n)); eb = np.array([st(rng.poisson(1.0, n).astype(float)) for _ in range(300)]).std(0, ddof=1)
    nul = np.array([((Lp := rng.permutation(Lm)) @ C / Lp.sum() - C.mean(0)) / BW for _ in range(200)]); en_ = nul.std(0, ddof=1); e = np.maximum(eb, en_)
    pr, vz, nv = NP.predikce(l, a); po = predikce_chybna(l, a)
    c2 = float(np.sum(((d - pr) / e) ** 2)); c2b = float(np.sum(((d - pr) / eb) ** 2)); a_, ea = amp(d - pr, e); c2o = float(np.sum(((d - po) / e) ** 2)); ao, _ = amp(d - po, e)
    out['skup'][nm] = dict(n=n, l=float(l.mean()), L_prumer=float(Lm.mean()), rozdil=d.tolist(), pred=pr.tolist(), pred_vazena=vz.tolist(), pred_nevazena=nv.tolist(),
                           vaz_dat=((Lm @ C / Lm.sum()) / BW).tolist(), nev_dat=(C.mean(0) / BW).tolist(), chyba_boot=eb.tolist(), chyba_nul=en_.tolist(), chyba=e.tolist(),
                           chi2=c2, p=float(stats.chi2.sf(c2, 5)), chi2_boot=c2b, nul_prumer=nul.mean(0).tolist(), a=a_, e_a=ea, mx1=x1,
                           pred_chybna=po.tolist(), chi2_chybna=c2o, p_chybna=float(stats.chi2.sf(c2o, 5)), a_chybna=ao)
    print('%-4s n %d l %.3f L̄ %.4f\n  rozdíl %s\n  pred   %s\n  chyba  %s\n  χ² %.2f (p %.3g; jen bootstrap %.2f)  amplituda %+.3f ± %.3f  ⟨x₁⟩ %.4f  |  chybná předpověď: χ² %.1f, amplituda %+.3f' % (nm, n, l.mean(), Lm.mean(), np.round(d, 5), np.round(pr, 5), np.round(e, 5), c2, out['skup'][nm]['p'], c2b, a_, ea, x1, c2o, ao), flush=True)
K, Z = out['skup']['even'], out['skup']['odd']
pmin = min(K['p'], Z['p']); out['verd_s7k'] = 'O' if (K['p'] > 1e-3 and Z['p'] > 1e-3) else ('Z' if pmin < 1e-4 else 'N')
# část B tak, jak je zapsaná, s opravenou předpovědí (doplněk D)
A0, EA0 = 0.608, 0.045; s_ = math.sqrt(K['e_a'] ** 2 + EA0 ** 2); zA = K['a'] / K['e_a']; dA = abs(K['a'] - A0)
out['verd_s7'] = 'O' if (zA > 5 and dA <= 2 * s_ and Z['p'] > 1e-3) else ('Z' if (dA > 3 * s_ or zA < 2) else 'N'); out['s7_z'] = zA; out['s7_odchylka_v_chybach'] = dA / s_
r = np.array(K['rozdil']) - np.array(K['pred']) - K['a'] * u; et = np.sqrt(np.array(K['chyba']) ** 2 + (K['a'] * eu) ** 2); c2t = float(np.sum((r / et) ** 2)); pt = float(stats.chi2.sf(c2t, 4))
out['tvar'] = dict(chi2=c2t, p=pt, verd='O' if pt > 0.01 else ('Z' if pt < 0.001 else 'N'), poznamka='amplituda je v chybě nula — test tvaru nemá co porovnávat (předpoklad padl)')
print('odv-lv-s7k %s (p sudé %.3g, liché %.3g)  |  odv-lv-s7 %s (a/chyba %.2f, |a − 0,608| = %.1f chyby)  |  odv-lv-tvar χ² %.2f p %.3g (%s, prázdné: a ≈ 0)' %
      (out['verd_s7k'], K['p'], Z['p'], out['verd_s7'], zA, dA / s_, c2t, pt, out['tvar']['verd']), flush=True)
json.dump(out, open(OUT, 'w'), indent=1); print('hotovo', round(time.time() - t0), 's')
