"""34.85a (popisně, bez verdiktu): proč je χ² vážené hustoty pod počtem košů? S7 rozdělená na K po sobě jdoucích kusů (podle |D|) a na K náhodných kusů;
rozptyl zbytku (data − opravená předpověď kusu) mezi kusy proti chybě, kterou pro kus dává bootstrap a promíchání. Kdyby byla chyba správná, je
⟨(zbytek/chyba)²⟩ ≈ 1 v každém koši. usage: python3 -I s7_kusy.py SLOZKA_NUL HODNOTY.npz VYSTUP.json [K]"""
import sys, os, gzip, json, math, numpy as np, scipy.special as sp
sys.path.insert(0, '/home/claude/br85/py'); import nez_porovnani as NP
ZD, HV, OUT = sys.argv[1:4]; K = int(sys.argv[4]) if len(sys.argv) > 4 else 10; rng = np.random.default_rng(34853); BW = 0.2
Dl, cnt, gl = [], [], []
for fn in sorted(os.listdir(ZD)):
    if fn.endswith('.jsonl.gz'):
        for line in gzip.open(os.path.join(ZD, fn), 'rt'):
            r = json.loads(line); Dl.append(int(r['D'])); g = r['g']; cnt.append(len(g)); gl.extend(g)
D = np.array(Dl, np.int64); cnt = np.array(cnt); G = np.array(gl, float); own = np.repeat(np.arange(len(D)), cnt)
z = np.load(HV); s = np.argsort(z['D']); HD, HL = z['D'][s], z['L'][s]; L = HL[np.searchsorted(HD, D)]
out = {'K': K}
for nm, sg, a in (('even', 1, 0), ('odd', -1, 1)):
    m = (D > 0) == (sg > 0); n = int(m.sum()); pos = np.full(len(D), -1); pos[m] = np.arange(n); Dm = np.abs(D[m]).astype(float); Lm = L[m]; l = np.log(Dm / math.pi)
    k = pos[own]; ok = k >= 0; t = G[ok]; k = k[ok]; x = ((t / 2) * l[k] + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi
    b = np.floor(x / BW).astype(int); inb = (x >= 0) & (b < 5); C = np.zeros((n, 5)); np.add.at(C, (k[inb], b[inb]), 1.0)
    st = lambda i, w=None: (((Lm[i] if w is None else w * Lm[i]) @ C[i]) / np.sum(Lm[i] if w is None else w * Lm[i]) - ((C[i].mean(0)) if w is None else (w @ C[i]) / w.sum())) / BW
    res = {}
    for jm, poradi in (('souvisle', np.argsort(Dm)), ('nahodne', rng.permutation(n))):
        ku = np.array_split(poradi, K); Z2 = []; ZB = []
        for i in ku:
            d = st(i); pr = NP.predikce(l[i], a)[0]
            eb = np.array([st(i, rng.poisson(1.0, len(i)).astype(float)) for _ in range(120)]).std(0, ddof=1)
            en = np.array([((Lp := rng.permutation(Lm[i])) @ C[i] / Lp.sum() - C[i].mean(0)) / BW for _ in range(80)]).std(0, ddof=1)
            Z2.append(((d - pr) / np.maximum(eb, en)) ** 2); ZB.append(((d - pr) / eb) ** 2)
        Z2 = np.array(Z2); ZB = np.array(ZB)
        res[jm] = dict(chi2_kusu=Z2.sum(1).tolist(), prumer_chi2=float(Z2.sum(1).mean()), prumer_chi2_jen_bootstrap=float(ZB.sum(1).mean()), prumer_z2_po_kosich=Z2.mean(0).tolist())
        print('%-4s %-9s K = %d: průměrné χ² kusu %.2f (jen bootstrap %.2f; čekáno 5 při správné chybě), po koších ⟨z²⟩ %s' % (nm, jm, K, res[jm]['prumer_chi2'], res[jm]['prumer_chi2_jen_bootstrap'], np.round(Z2.mean(0), 2)), flush=True)
    out[nm] = res
json.dump(out, open(OUT, 'w'), indent=1)
