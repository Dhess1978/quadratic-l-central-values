"""34.79 doplněk k odv-lr-vaha (David 12. 11. se ptal, zda odchylka není chybou přerušeného kontrolního balíku).
(1) data rozdílu vážené − nevážené hustoty znovu přímo ze surových balíků nul (KS jsonl.gz) a hodnot L (V/D100), nezávisle na h3_hust.py;
(2) nulový model podle pravidla: promíchání L mezi D (200×) → rozdíl ≈ 0, χ² proti 0 a proti předpovědi; rozptyl promíchání proti chybám bootstrapu;
(3) rozklad sudé/liché podle třídy D: 2 ∤ D (χ(2) = +1, D ≡ 1 mod 8; χ(2) = −1, D ≡ 5 mod 8) a 4 | D."""
import os, gzip, json, math, numpy as np, scipy.special as sp
from scipy import stats
SRC = {'B10': '/home/claude/br63/KS/laborator_male_2026_11_02/B10', 'C100': '/home/claude/br63/KS/velky_2026_11_02/C100', 'R10': '/home/claude/br63/KS/velky_2026_11_02/R10', 'D100': '/home/claude/br63/KS/velky_2026_11_02/D100'}
FZ = json.load(open('/home/claude/br79/data/fazzari.json')); H3 = json.load(open('/home/claude/br78/data/h3_hust.json'))
EDG = np.linspace(0, 1, 6); BW = 0.2; rng = np.random.default_rng(34791); out = {}
for F, P in SRC.items():
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
    Ds, Gs = [], []
    for fn in sorted(os.listdir(P)):
        if fn.endswith('.jsonl.gz'):
            for line in gzip.open(os.path.join(P, fn), 'rt'): r = json.loads(line); Ds.append(int(r['D'])); Gs.append(r['g'])
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        idx = [i for i, d in enumerate(Ds) if (d > 0) == (sg > 0) and d in lv]; n = len(idx)
        D = np.array([Ds[i] for i in idx], np.int64); L = np.array([lv[Ds[i]] for i in idx]); C = np.zeros((n, 5)); n0 = 0
        for j, i in enumerate(idx):
            t = np.asarray(Gs[i], float)
            if t.size == 0: n0 += 1; continue
            x = ((t / 2) * math.log(abs(Ds[i]) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi; C[j] = np.histogram(x, EDG)[0]
        roz = lambda w: ((w * L) @ C / np.sum(w * L) - w @ C / np.sum(w)) / BW
        d = roz(np.ones(n)); key = F + '_' + nm; g = FZ['skupiny'][key]; e = np.array(g['chyba']); pred = np.array(g['pred'])
        nul = np.array([((Lp := rng.permutation(L)) @ C / Lp.sum() - C.sum(0) / n) / BW for _ in range(200)]); sn = nul.std(0, ddof=1)
        c2n = [float(np.sum((v / e) ** 2)) for v in nul]
        # třídy D
        tr = {}
        for tn, m in (('1mod8', (np.abs(D) % 8 == 1) if sg > 0 else (D % 8 == 1)), ('5mod8', (D % 8 == 5)), ('4delí', (D % 4 == 0))):
            if m.sum() < 100: continue
            Cm, Lm = C[m], L[m]; dm = ((Lm @ Cm) / Lm.sum() - Cm.sum(0) / m.sum()) / BW
            B = np.array([(lambda w: ((w * Lm) @ Cm / np.sum(w * Lm) - w @ Cm / np.sum(w)) / BW)(rng.poisson(1.0, m.sum()).astype(float)) for _ in range(150)]); em = B.std(0, ddof=1)
            tr[tn] = dict(n=int(m.sum()), rozdil=dm.tolist(), chyba=em.tolist(), chi2_pred=float(np.sum(((dm - pred) / em) ** 2)), p=float(stats.chi2.sf(float(np.sum(((dm - pred) / em) ** 2)), 5)))
        out[key] = dict(n=n, bez_nul=n0, rozdil_znovu=d.tolist(), rozdil_h3=H3['skupiny'][key]['rozdil'][:5], max_rozdil_h3=float(np.max(np.abs(d - np.array(H3['skupiny'][key]['rozdil'][:5])))),
                        nul_prumer=nul.mean(0).tolist(), nul_sd=sn.tolist(), chyba_boot=e.tolist(), pomer_sd=(sn / e).tolist(), chi2_nul_median=float(np.median(c2n)), chi2_nul_95=float(np.quantile(c2n, 0.95)),
                        chi2_nul_proti_pred=float(np.sum(((nul.mean(0) - pred) / e) ** 2)), tridy=tr)
        print('%-9s n %6d (bez nul %d)  znovu − h3 max %.1e | nul: průměr %s sd/boot %s χ² medián %.1f (95 %% %.1f) | třídy %s' % (key, n, n0, out[key]['max_rozdil_h3'], np.round(nul.mean(0), 4), np.round(sn / e, 2), np.median(c2n), np.quantile(c2n, 0.95),
              {k: (v['n'], round(v['chi2_pred'], 1)) for k, v in tr.items()}), flush=True)
json.dump(out, open('/home/claude/br79/data/vaha_nul.json', 'w'), indent=1)
