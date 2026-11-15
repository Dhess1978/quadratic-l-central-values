"""34.82: L(½, χ_D) vlastní vyhlazenou funkcionální rovnicí (numpy), obnovitelně po dávkách; váhy Q(c, πn²/|D|) Taylorem 2. řádu kolem
|D|₀ = 1,0001·10¹⁰ (|D| se v E10 + E100 liší o méně než 10⁻⁵ relativně; zbytek ~ (55·10⁻⁵)³ — zanedbatelný).
L(½, χ_D) = 2 Σ_{n ≥ 1} χ_D(n) n^{−½} Q(c, πn²/|D|), Q = horní regularizovaná neúplná Γ, c = ¼ (D > 0) a ¾ (D < 0), n ≤ √(55|D|/π) + 2.
χ_D(n) úplně multiplikativně z χ_D(p): p = 2 podle D mod 8, liché p Eulerovým kritériem; tabulka nejmenších prvočinitelů.
Použití: python3 -I afe_L.py vstup.npz vystup_prefix cast pocet_casti"""
import sys, os, math, time, numpy as np
from scipy.special import gammaincc
vst, pref, cast, casti = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
D_all = np.load(vst)['D'].astype(np.int64); idx = np.arange(len(D_all))[cast::casti]
NMAX = int(math.sqrt(55 * int(np.abs(D_all).max()) / math.pi)) + 3
spf = np.zeros(NMAX + 1, np.int64)
for p in range(2, int(NMAX ** 0.5) + 1):
    if spf[p] == 0: m = spf[p * p::p] == 0; spf[p * p::p][m] = p
n = np.arange(NMAX + 1); prime = (spf == 0) & (n >= 2); spf[prime] = n[prime]
P = np.nonzero(prime)[0]; pidx = np.full(NMAX + 1, -1, np.int64); pidx[P] = np.arange(len(P))
nn = np.arange(2, NMAX + 1); sp = spf[2:]; q = nn // sp; ip = pidx[sp]
E = (P - 1) // 2; nb = int(E.max()).bit_length()
from scipy.special import gamma as _G
Q0 = 1.0001e10; kk = np.arange(1, NMAX + 1).astype(float); x0 = math.pi * kk ** 2 / Q0; TAB = {}
for c in (0.25, 0.75):
    g = np.exp((c - 1) * np.log(x0) - x0) / _G(c); TAB[c] = (gammaincc(c, x0), g, g * ((c - 1) / x0 - 1), x0, 1 / np.sqrt(kk))   # Q, −dQ/dx, d(−dQ/dx)/dx
def chi_p(D):
    r = D % P; res = np.ones(len(P), np.int64); b = r.copy(); e = E.copy()
    for _ in range(nb):
        m = (e & 1) == 1; res[m] = (res[m] * b[m]) % P[m]; b = (b * b) % P; e >>= 1
    c = np.where(r == 0, 0, np.where(res == 1, 1, -1)); c[0] = 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1); return c.astype(np.int8)
OUT = pref + '_%d_%d.npz' % (cast, casti); hot = dict(np.load(OUT)) if os.path.exists(OUT) else {'i': np.zeros(0, np.int64), 'L': np.zeros(0)}
done = set(hot['i'].tolist()); I, Lv = list(hot['i']), list(hot['L']); t0 = time.time(); k0 = len(I)
for j in idx:
    if j in done: continue
    D = int(D_all[j]); aq = abs(D); N = int(math.sqrt(55 * aq / math.pi)) + 2; c = 0.25 if D > 0 else 0.75
    cp = chi_p(D); ch = np.zeros(N + 1, np.int8); ch[1] = 1; a = ip[:N - 1]; qq = q[:N - 1]
    for _ in range(25):
        new = cp[a] * ch[qq]
        if np.array_equal(new, ch[2:N + 1]): break
        ch[2:N + 1] = new
    T = TAB[c]; h = -T[3][:N] * (aq - Q0) / aq          # x = x0·Q0/q = x0(1 − u), h = −x0·u
    w = (T[0][:N] - T[1][:N] * h - 0.5 * T[2][:N] * h * h) * T[4][:N]
    I.append(j); Lv.append(2.0 * float(np.dot(ch[1:].astype(float), w)))
    if (len(I) - k0) % 1000 == 0:
        np.savez(OUT, i=np.array(I), L=np.array(Lv)); print(cast, len(I), round(time.time() - t0), 's', flush=True)
np.savez(OUT, i=np.array(I), L=np.array(Lv)); print('hotovo', cast, len(I), round(time.time() - t0), 's', flush=True)
