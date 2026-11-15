"""34.85 odv-lv-ef a odv-lv-hladka: explicitní vzorec na celé rodině (strana nul proti straně prvočísel z naměřených L(½)) a srovnání s předpovědí.
φ(u) = (sin(πu/2)/(πu/2))^4, u = γ ln(|D|/π)/(2π), φ̂ nosič |ν| ≤ 1 → n ≤ |D|/π. Konvence ověřeny na jednotlivých L (py/ef_jeden.py, shoda 1e-6).
usage: python3 ef.py RODINA ZNAMENKO(even|odd) ADRESAR_NUL[,ADRESAR2] HODNOTY.npz[,…] VYSTUP.json"""
import sys, os, gzip, json, math, time, numpy as np, scipy.special as sp
from scipy import integrate
sys.path.insert(0, '/home/claude/br85/py')
t0 = time.time(); F, ZN, ZD, HV, OUT = sys.argv[1:6]; sg, a = (1, 0) if ZN == 'even' else (-1, 1)
phi = lambda u: np.sinc(u / 2) ** 4
def M4(x):
    x = np.abs(x); return np.where(x <= 1, 2 / 3 - x ** 2 + x ** 3 / 2, np.where(x <= 2, (2 - x) ** 3 / 6, 0.0))
phihat = lambda nu: 2 * M4(2 * nu)
# ---- nuly
Ds, Zp, Gm = [], [], []
for zd in ZD.split(','):
    for fn in sorted(os.listdir(zd)):
        if not fn.endswith('.jsonl.gz'): continue
        for line in gzip.open(os.path.join(zd, fn), 'rt'):
            r = json.loads(line); D = int(r['D'])
            if (D > 0) != (sg > 0): continue
            g = np.asarray(r['g'], float); l = math.log(abs(D) / math.pi)
            Ds.append(D); Zp.append(phi(g * l / (2 * math.pi)).sum()); Gm.append(g.max() if len(g) else 0.0)
o = np.argsort(np.abs(Ds)); D = np.array(Ds, np.int64)[o]; Zp = np.array(Zp)[o]; N = len(D); print(F, ZN, 'L-funkcí', N, round(time.time() - t0), 's', flush=True)
# ---- hodnoty
HD, HL = [], []
for f in HV.split(','):
    z = np.load(f); HD.append(z['D']); HL.append(z['L'])
HD = np.concatenate(HD); HL = np.concatenate(HL); s = np.argsort(HD); HD, HL = HD[s], HL[s]; ix = np.searchsorted(HD, D); assert np.all(HD[ix] == D), 'hodnoty chybí'; L = HL[ix]
l = np.log(np.abs(D) / math.pi)
# ---- archimedovský člen a ocas nad spočtenými nulami (x > 4,35) na mřížce l
lg = np.linspace(l.min() - 1e-6, l.max() + 1e-6, 41)
def arch(lv): return 4 / 3 + integrate.quad(lambda u: phi(u) * np.real(sp.digamma((0.5 + a + 2j * math.pi * u / lv) / 2)), -400, 400, limit=8000, points=[0])[0] / lv
def theta(t, lv): return 0.5 * t * lv + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))
def tmax(lv, x=4.35):
    T = 2.0
    for _ in range(60):
        f = theta(T, lv) / math.pi - x; d = (0.5 * lv + 0.5 * np.real(sp.digamma((0.5 + a + 1j * T) / 2))) / math.pi; T -= f / d
    return T
def ocas(lv):
    T = tmax(lv); f = lambda t: phi(t * lv / (2 * math.pi)) * (0.5 * lv + 0.5 * np.real(sp.digamma((0.5 + a + 1j * t) / 2))) / math.pi
    return integrate.quad(f, T, T + 600, limit=4000)[0]
AR = np.interp(l, lg, [arch(v) for v in lg]); OC = np.interp(l, lg, [ocas(v) for v in lg])
Z = 2 * (Zp + OC)                     # celá přímka (nuly ±γ)
print('arch, ocas', round(time.time() - t0), 's; ocas průměr %.2e' % OC.mean(), flush=True)
# ---- prvočíselná strana: −(2/l) Σ Λ(n) χ_D(n) n^{-1/2} φ̂(ln n / l), rozložená na části
nmax = int(np.abs(D).max() / math.pi); sv = np.ones(nmax + 1, bool); sv[:2] = False
for q in range(2, int(nmax ** 0.5) + 1):
    if sv[q]: sv[q * q::q] = False
PR = np.nonzero(sv)[0]
parts = {k: np.zeros(N) for k in ('p2', 'kratke', 'dlouhe', 'mocniny')}
Dm8 = D % 8; chi2 = np.where(D % 2 == 0, 0, np.where((Dm8 == 1) | (Dm8 == 7), 1, -1)).astype(float)
for j, p in enumerate(PR):
    lp = math.log(p)
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
    if j % 5000 == 0: print('prvočísla', j, len(PR), round(time.time() - t0), 's', flush=True)
Pp = sum(parts.values()); P = AR + Pp
res = P - Z
print('rozdíl na jednotlivých L: medián %.2e, 99,9 %% |·| %.2e, max |·| %.2e' % (np.median(res), np.quantile(np.abs(res), 0.999), np.abs(res).max()), flush=True)
# ---- vážené − nevážené a bootstrap
rng = np.random.default_rng(34850); B = 300
def delta(X, w): return (w * L) @ X / (w @ L) - w @ X / w.sum()
Q = {'Z': Z, 'P': P, 'arch': AR, 'ocas': 2 * OC, **{'P_' + k: v for k, v in parts.items()}}
est = {k: float(delta(v, np.ones(N))) for k, v in Q.items()}; bs = {k: [] for k in Q}; bs['Z-P'] = []
for b in range(B):
    w = rng.poisson(1.0, N).astype(float)
    for k, v in Q.items(): bs[k].append(delta(v, w))
    bs['Z-P'].append(bs['Z'][-1] - bs['P'][-1])
err = {k: float(np.std(v, ddof=1)) for k, v in bs.items()}; est['Z-P'] = est['Z'] - est['P']
out = dict(rodina=F, znamenko=ZN, n=int(N), l=float(l.mean()), delta=est, chyba=err, rozdil_L=dict(median=float(np.median(res)), q999=float(np.quantile(np.abs(res), 0.999)), max=float(np.abs(res).max()),
           nejhorsi=[[int(D[i]), float(res[i])] for i in np.argsort(-np.abs(res))[:10]]), M=float(L.mean()))
# ---- předpověď (nezávislá implementace) pro tutéž φ
import nez_porovnani as NP
TG, DT = NP.TG, NP.DT
ed = np.linspace(l.min(), l.max() + 1e-9, 41); idx = np.digitize(l, ed) - 1; nw = nm = vw = vu = 0.0
for i in range(40):
    m = idx == i
    if not m.any(): continue
    w = m.sum(); lc = l[m].mean(); Mv, r0, s_w, s_u = NP.hustoty(lc, a); ph = phi(TG * lc / (2 * math.pi))
    vw += w * np.sum(ph * (Mv * r0 + s_w)) * DT; nm += w * Mv; vu += w * np.sum(ph * (r0 + s_u)) * DT; nw += w
out['delta_pred'] = 2 * (vw / nm - vu / nw)
out['hladka_z'] = (est['Z'] - out['delta_pred']) / err['Z']; out['ef_z'] = est['Z-P'] / err['Z-P']
print(json.dumps({k: out[k] for k in ('n', 'l', 'delta', 'chyba', 'delta_pred', 'hladka_z', 'ef_z')}, indent=None)[:1500])
np.savez(OUT.replace('.json', '_D.npz'), D=D, L=L, Z=Z, P=P, OC=OC, **{'P_' + k: v for k, v in parts.items()})
for mez in (1e-2, 3e-3):   # popisně: bez L-funkcí s velkým rozdílem na jednotlivé L (podezření na vynechanou blízkou dvojici nul, divz 8)
    keep = np.abs(res) < mez; zz = []; rng2 = np.random.default_rng(34851)
    for b in range(100):
        w = rng2.poisson(1.0, N).astype(float) * keep; zz.append(delta(Z - P, w))
    out['bez_velkych_%g' % mez] = dict(vynechano=int((~keep).sum()), Z_P=float(delta(Z - P, keep.astype(float))), chyba=float(np.std(zz, ddof=1)))
    print('bez |rozdíl| >= %g (%d L): Z − P = %.2e ± %.1e' % (mez, (~keep).sum(), out['bez_velkych_%g' % mez]['Z_P'], out['bez_velkych_%g' % mez]['chyba']), flush=True)
json.dump(out, open(OUT, 'w'), indent=1); print('hotovo', round(time.time() - t0), 's')
