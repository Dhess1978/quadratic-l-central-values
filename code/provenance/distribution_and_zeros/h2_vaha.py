"""34.78 H2 odv-lq-vaha: posun střední první nuly při vážení L(½)^k proti USp(2N_x1) vážené Λ(1)^k (přesně, Jacobiho soubor)."""
import sys, csv, math, json, numpy as np
sys.path.insert(0, '/home/claude/br78/py'); sys.path.insert(0, '/home/claude/br76/py'); import usp_jacobi as UJ, usp_dpp
# 0) nezávislá kontrola přesného modelu: DPP vzorky převážené Λ
lnL, x1 = usp_dpp.sample(3, 20000, seed=378); w1 = np.exp(lnL); w2 = w1 ** 2
kon = dict(N=3, dpp_k1=float(np.sum(w1 * x1) / np.sum(w1)), dpp_k2=float(np.sum(w2 * x1) / np.sum(w2)), presne_k1=UJ.ex1(3, 1)[0], presne_k2=UJ.ex1(3, 2)[0])
print('kontrola DPP (20 000) proti přesnému:', kon, flush=True)
# 1) přesné USp: E[x1] pro k = 0, 1, 2 a N = 1..12
U = {N: [UJ.ex1(N, k)[0] for k in (0, 1, 2)] for N in range(1, 13)}
def interp(m0):  # N_x1 z neváženého E[x1] (lineárně mezi celými N), vrátí N a Δ_1, Δ_2
    Ns = sorted(U); m = [U[N][0] for N in Ns]; i = np.searchsorted(m, m0) - 1; i = min(max(i, 0), len(Ns) - 2); t = (m0 - m[i]) / (m[i + 1] - m[i])
    d = [(U[Ns[i]][k] - U[Ns[i]][0]) * (1 - t) + (U[Ns[i + 1]][k] - U[Ns[i + 1]][0]) * t for k in (1, 2)]
    return Ns[i] + t, d
# 2) data: první nuly (surové γ) a L, rozvinutí θ_D jako ve studii III / 34.77
def lngamma_c(z):
    s = np.zeros_like(z)
    for k in range(12): s = s - np.log(z + k)
    w = z + 12; B = [1/12, -1/360, 1/1260, -1/1680, 1/1188, -691/360360, 1/156, -3617/122400]
    t = (w - 0.5) * np.log(w) - w + 0.5 * math.log(2 * math.pi); iw = 1 / w; p = iw
    for b in B: t = t + b * p; p = p * iw * iw
    return s + t
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as fh:
    rd = csv.reader(fh); next(rd)
    for F, d0, g1, g2, nn in rd: rows.setdefault(F, []).append((int(d0), float(g1)))
rng = np.random.default_rng(3478); out = dict(kontrola_dpp=kon, usp={N: U[N] for N in U}, skupiny={})
for F in ('B10', 'C100', 'R10', 'D100'):
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
    for nm, s, a in (('kl', 1, 0), ('zap', -1, 1)):
        pr = [(d0, g1, lv[d0]) for d0, g1 in rows[F] if (d0 > 0) == (s > 0) and d0 in lv and math.isfinite(g1)]
        Dv = np.array([p[0] for p in pr], np.int64); t = np.array([p[1] for p in pr]); L = np.array([p[2] for p in pr])
        x = ((t / 2) * np.log(np.abs(Dv) / math.pi) + np.imag(lngamma_c((0.5 + a + 1j * t) / 2))) / math.pi
        def stat(L, x): m0 = x.mean(); return [m0] + [float(np.sum(L ** k * x) / np.sum(L ** k) - m0) for k in (1, 2)]
        S = stat(L, x); n = len(x); B = np.array([stat(L[i], x[i]) for i in (rng.integers(0, n, n) for _ in range(1000))])
        Nx, dU = interp(S[0]); e = B.std(axis=0, ddof=1)
        sh = [stat(rng.permutation(L), x) for _ in range(20)]
        g = dict(n=n, mx1=S[0], N_x1=Nx, d=S[1:], e=e[1:].tolist(), d_usp=dU, z=[(S[1 + k] - dU[k]) / e[1 + k] for k in (0, 1)], nul_d=np.mean([q[1:] for q in sh], axis=0).tolist())
        out['skupiny']['%s_%s' % (F, nm)] = g
        print('%-9s n %6d  E[x1] %.4f  N_x1 %.2f  Δ1 %.4f ± %.4f (USp %.4f, z %+.1f)  Δ2 %.4f ± %.4f (USp %.4f, z %+.1f)  nulový %.4f %.4f' % ('%s_%s' % (F, nm), n, S[0], Nx, S[1], e[1], dU[0], g['z'][0], S[2], e[2], dU[1], g['z'][1], *g['nul_d']), flush=True)
zmax = max(abs(z) for g in out['skupiny'].values() for z in g['z']); out['verdikt'] = 'O' if zmax < 3 else ('Z' if zmax >= 4 else 'N'); out['zmax'] = zmax
print('max |z|', round(zmax, 2), 'verdikt', out['verdikt'])
json.dump(out, open('/home/claude/br78/data/h2_vaha.json', 'w'), indent=1)
