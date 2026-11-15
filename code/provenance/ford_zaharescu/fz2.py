"""34.84 · Ford–Zaharescu: rozdělení {α·γₙ} mod 1 podle PRAVIDLA_3484.md (oddíl 3, zapsáno před výpočtem).
Předpověď bez volného parametru z Landauova vzorce: f(u) = 1 − 2w·Σ_k Λ(p^k)·p^{−k/2}·cos(2πku), w = (T₂ − T₁)/(2πN).
Data: 100 000 nul LMFDB (1 blok), zeros6 ve 4 blocích, Odlyzkovy bloky 10¹², 10²¹, 10²² (popisně).
Nulový model: tytéž rozestupy v náhodném pořadí (200×). Chyba koše = max(√očekávání, sd nulového modelu)."""
import json, math, re, numpy as np, mpmath as mp
from scipy import stats
rng = np.random.default_rng(34841); NB = 50; ED = np.linspace(0, 1, NB + 1); KMAX = 60; NNUL = 200
mp.mp.dps = 40
def nacti_lmfdb():
    g = np.array([float(l.split()[1]) for l in open('/home/claude/lzin14/zetazeros100000.txt') if l.strip()]); return g
g6 = np.loadtxt('/home/claude/z6/zeros6.txt')
BL = [('L1e5', nacti_lmfdb(), None)]
q = len(g6) // 4
for i in range(4): BL.append(('Z6%s' % 'abcd'[i], g6[i * q:(i + 1) * q] if i < 3 else g6[3 * q:], None))
for nm, f in (('O1e12', 'zeros3'), ('O1e21', 'zeros4'), ('O1e22', 'zeros5')):
    T = open('/home/claude/pub4/run1/odl/%s.txt' % f).read().split('\n'); B = int(re.search(r'gamma - (\d+)', T[0]).group(1))
    d = np.array([float(l) for l in T if re.fullmatch(r'\s*\d+\.\d+\s*', l)]); BL.append((nm, d, B))
ALFY = [('p2', 2, True), ('p3', 3, True), ('p5', 5, True), ('p7', 7, True), ('odm2', None, False), ('ln6', 6, False), ('ln10', 10, False), ('ln4', 4, 'moc')]
def alfa(a):
    if a[0] == 'odm2': return mp.sqrt(2) / (2 * mp.pi)
    return mp.log(a[1]) / (2 * mp.pi)
def zlomky(g, B, al):
    if B is None: return np.mod(float(al) * g, 1.0)
    f0 = float(mp.frac(al * B)); return np.mod(f0 + float(al) * g, 1.0)
def pred_kose(w, a):
    """očekávaný podíl v koších: (b − a) − 2w Σ_k c_k (sin 2πkb − sin 2πka)/(2πk)"""
    c = np.zeros(KMAX + 1)
    if a[2] is True: c[1:] = [math.log(a[1]) / a[1] ** (k / 2) for k in range(1, KMAX + 1)]
    elif a[2] == 'moc': c[1:] = [math.log(2) / 2 ** k for k in range(1, KMAX + 1)]   # 4^k = 2^{2k}: Λ = ln 2, (4^k)^{1/2} = 2^k
    P = np.diff(ED).copy()
    for k in range(1, KMAX + 1):
        if c[k]: P -= 2 * w * c[k] * (np.sin(2 * np.pi * k * ED[1:]) - np.sin(2 * np.pi * k * ED[:-1])) / (2 * np.pi * k)
    return P, c
out = {'bloky': {}, 'prav': {}}
for nm, g, B in BL:
    N = len(g); sp = np.diff(g); h = sp.mean(); T1 = (g[0] - h / 2) if not (nm == 'L1e5') else 0.0; T2 = g[-1] + h / 2
    if B is not None: T1 += 0; Tm = B + (g[0] + g[-1]) / 2
    else: Tm = (g[0] + g[-1]) / 2
    w = (T2 - T1) / (2 * math.pi * N); out['bloky'][nm] = dict(N=N, w=w, t1=float(g[0]) if B is None else None, t2=float(g[-1]) if B is None else None, zaklad=B, lnT=float(math.log(Tm / (2 * math.pi))), alfy={})
    for a in ALFY:
        al = alfa(a); u = zlomky(g, B, al); O = np.histogram(u, ED)[0].astype(float); P, c = pred_kose(w, a); E = N * P
        Nh = np.zeros((NNUL, NB))
        for r in range(NNUL):
            gg = g[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(sp))]); Nh[r] = np.histogram(zlomky(gg, B, al), ED)[0]
        sdn = Nh.std(0, ddof=1); err = np.maximum(np.sqrt(E), sdn); Eu = np.full(NB, N / NB); erru = np.maximum(np.sqrt(Eu), sdn)
        c2p = float(np.sum(((O - E) / err) ** 2)); c2u = float(np.sum(((O - Eu) / erru) ** 2))
        c1 = float(np.mean(np.cos(2 * np.pi * u))); s1 = float(np.mean(np.sin(2 * np.pi * u))); c1n = Nh @ np.cos(2 * np.pi * (ED[:-1] + ED[1:]) / 2) / N
        out['bloky'][nm]['alfy'][a[0]] = dict(chi2_pred=c2p, p_pred=float(stats.chi2.sf(c2p, NB - 1)), chi2_rovn=c2u, p_rovn=float(stats.chi2.sf(c2u, NB - 1)),
            c1=c1, s1=s1, c1_pred=float(-w * c[1]), c1_sd_nul=float(c1n.std(ddof=1)), sd_nul_rel=float(np.mean(sdn / np.sqrt(Eu))), O=O.tolist(), E=E.tolist(), err=err.tolist(),
            amp1=float(2 * w * c[1]), sdn=sdn.tolist(), c1n=c1n.tolist())
        print('%-6s %-5s χ²/49 předpověď %7.1f (p %.3g)  rovnoměrné %9.1f (p %.3g)  c1 %+.5f pred %+.5f (sd nul %.5f)' % (nm, a[0], c2p, out['bloky'][nm]['alfy'][a[0]]['p_pred'], c2u, out['bloky'][nm]['alfy'][a[0]]['p_rovn'], c1, -w * c[1], c1n.std(ddof=1)), flush=True)
HL = ['L1e5', 'Z6a', 'Z6b', 'Z6c', 'Z6d']
pp = [out['bloky'][b]['alfy'][a]['p_pred'] for b in HL for a in ('p2', 'p3', 'p5', 'p7')]
pu2 = [out['bloky'][b]['alfy']['p2']['p_rovn'] for b in HL[1:]]
out['prav']['odv-fz-tvar'] = 'O' if (min(pp) > 1e-3 and max(pu2) < 1e-6) else ('Z' if min(pp) < 1e-4 else 'N')
pk = [out['bloky'][b]['alfy'][a]['p_rovn'] for b in HL for a in ('odm2', 'ln6', 'ln10')]
out['prav']['odv-fz-kontrola'] = 'O' if min(pk) > 1e-3 else ('Z' if min(pk) < 1e-4 else 'N')
out['min_p_tvar'] = min(pp); out['max_p_rovn_p2'] = max(pu2); out['min_p_kontrola'] = min(pk)
print('verdikty', out['prav'], 'min p tvar %.3g, max p rovn (p=2) %.3g, min p kontrola %.3g' % (min(pp), max(pu2), min(pk)))
json.dump(out, open('/home/claude/br84/data/fz2.json', 'w'))
