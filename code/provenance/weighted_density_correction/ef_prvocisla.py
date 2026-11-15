"""34.85a: prvočíselná strana explicitního vzorce pro seznam D — stejné členy jako ef.py (−(2/l) Σ Λ(n) χ_D(n) n^{−1/2} φ̂(ln n / l), části p = 2, krátká, dlouhá, mocniny),
jen rychlejší cesta pro velká prvočísla: pro p > √(max|D|/π) je k = 1 a část „dlouhá“, ν = ln p / l ∈ (½, 1), tedy φ̂(ν) = 2·M4(2ν) = (8/3)(1 − ν)³;
tabulka Legendreova symbolu v int8. Malá prvočísla (p ≤ P0) jdou původní smyčkou ef.py beze změny. Kontrola proti uloženým částem C100: --kontrola.
usage: python3 -I ef_prvocisla.py RODINA ZNAMENKO LO HI HODNOTY.npz[,…] ZASOBA.npz
       python3 -I ef_prvocisla.py --kontrola ZNAMENKO ef2_C100_ZN_D.npz [kolik D]"""
import sys, os, math, time, numpy as np
t0 = time.time()
def M4(x):
    x = np.abs(x); return np.where(x <= 1, 2 / 3 - x ** 2 + x ** 3 / 2, np.where(x <= 2, (2 - x) ** 3 / 6, 0.0))
phihat = lambda nu: 2 * M4(2 * nu)
def casti(D, tisk=True):
    """části strany prvočísel pro pole D (int64, jedno znaménko)"""
    N = len(D); l = np.log(np.abs(D) / math.pi); nmax = int(np.abs(D).max() / math.pi); sv = np.ones(nmax + 1, bool); sv[:2] = False
    for q in range(2, int(nmax ** 0.5) + 1):
        if sv[q]: sv[q * q::q] = False
    PR = np.nonzero(sv)[0]; P0 = int(nmax ** 0.5) + 2      # nad P0: p² > nmax (jen k = 1) a ln p > l/2 pro všechna D (jen „dlouhá“)
    parts = {k: np.zeros(N) for k in ('p2', 'kratke', 'dlouhe', 'mocniny')}
    Dm8 = D % 8; chi2 = np.where(D % 2 == 0, 0, np.where((Dm8 == 1) | (Dm8 == 7), 1, -1)).astype(float)
    invl = 1.0 / l; acc = np.zeros(N); w = np.empty(N); nP = len(PR)
    for j, p in enumerate(PR):
        p = int(p); lp = math.log(p)
        if p <= P0:      # původní smyčka ef.py
            if p == 2: ch = chi2
            else:
                tab = -np.ones(p); tab[(np.arange(1, (p - 1) // 2 + 1) ** 2) % p] = 1; tab[0] = 0; ch = tab[D % p]
            k = 1; pk = p
            while pk <= nmax:
                nu = k * lp / l; m = nu < 1
                if not m.any(): break
                c = ch if k % 2 else ch * ch
                term = np.zeros(N); term[m] = -(2 / l[m]) * lp * c[m] / math.sqrt(pk) * phihat(nu[m])
                if p == 2: parts['p2'] += term
                elif k >= 2: parts['mocniny'] += term
                else:
                    kr = lp <= l / 2; parts['kratke'] += np.where(kr, term, 0); parts['dlouhe'] += np.where(kr, 0, term)
                k += 1; pk *= p
        else:            # rychlá cesta: k = 1, dlouhá část, φ̂(ν) = (8/3)(1 − ν)³ pro ½ < ν < 1, jinak 0
            h = (p - 1) // 2; kk = np.arange(1, h + 1, dtype=np.int64); tab = np.full(p, -1, np.int8); tab[(kk * kk) % p] = 1; tab[0] = 0
            np.multiply(invl, -lp, out=w); w += 1.0; np.maximum(w, 0.0, out=w)      # (1 − ν)₊
            acc += tab[D % p] * (w * w * w) * (invl * (-(16.0 / 3.0) * lp / math.sqrt(p)))
        if tisk and j % 20000 == 0: print('prvočísla', j, nP, round(time.time() - t0), 's', flush=True)
    parts['dlouhe'] += acc; return parts
if sys.argv[1] == '--kontrola':
    ZN, FN = sys.argv[2:4]; K = int(sys.argv[4]) if len(sys.argv) > 4 else 20000; z = np.load(FN); rng = np.random.default_rng(1); i = np.sort(rng.choice(len(z['D']), K, replace=False)); i[-1] = len(z['D']) - 1
    D = z['D'][i].astype(np.int64); P = casti(D, False)     # poslední D = největší |D| → stejné nmax jako v ef.py
    for k in P: d = np.abs(P[k] - z['P_' + k][i]); print('%-8s největší rozdíl %.2e (hodnoty do %.3f)' % (k, d.max(), np.abs(z['P_' + k][i]).max()))
    print('hotovo', round(time.time() - t0), 's'); sys.exit(0)
F, ZN = sys.argv[1:3]; LO, HI = int(sys.argv[3]), int(sys.argv[4]); HV, ZAS = sys.argv[5:7]; sg = 1 if ZN == 'even' else -1
HD = np.concatenate([np.load(f)['D'] for f in HV.split(',')]); D = HD[(np.abs(HD) > LO) & (np.abs(HD) <= HI) & ((HD > 0) == (sg > 0))]; D = D[np.argsort(np.abs(D))].astype(np.int64)
print(F, ZN, 'D z tabulky hodnot', len(D), flush=True); P = casti(D)
np.savez(ZAS, D=D, **{'P_' + k: v for k, v in P.items()}); print('zásoba uložena', ZAS, len(D), 'D,', round(time.time() - t0), 's')
