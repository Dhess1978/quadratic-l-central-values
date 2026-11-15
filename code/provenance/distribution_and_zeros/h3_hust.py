"""34.78 H3 odv-lq-hust: rozdíl hustoty nízkých nul vážené L(½) a nevážené proti USp(2N_x1) vážené Λ(1) (přesně, Jacobiho soubor)."""
import sys, os, gzip, json, math, numpy as np
from scipy import stats
sys.path.insert(0, '/home/claude/br78/py'); import usp_jacobi as UJ
from h2_vaha_lng import lngamma_c
H2 = json.load(open('/home/claude/br78/data/h2_vaha.json'))
edges = np.linspace(0, 3, 16); nb = 15; bw = 0.2
def usp_rozdil(N):  # hustota v x = θN/π po koších: vážená (k=1) − nevážená
    th = np.linspace(1e-9, np.pi - 1e-9, 200001); out = []
    for k in (0, 1):
        h = UJ.hustota_theta(N, k, th) * np.pi / N; xx = th * N / np.pi
        c = np.concatenate([[0], np.cumsum((h[1:] + h[:-1]) / 2 * np.diff(xx))]); out.append(np.diff(np.interp(edges, xx, c)) / bw)
    return out[1] - out[0], out
R = {N: usp_rozdil(N) for N in range(1, 7)}
def usp_N(Nx): i = int(math.floor(Nx)); t = Nx - i; return R[i][0] * (1 - t) + R[i + 1][0] * t, [R[i][1][k] * (1 - t) + R[i + 1][1][k] * t for k in (0, 1)]
SRC = {'B10': '/home/claude/br63/KS/laborator_male_2026_11_02/B10', 'C100': '/home/claude/br63/KS/velky_2026_11_02/C100', 'R10': '/home/claude/br63/KS/velky_2026_11_02/R10', 'D100': '/home/claude/br63/KS/velky_2026_11_02/D100'}
rng = np.random.default_rng(34781); out = {'edges': edges.tolist(), 'skupiny': {}}
for F, P in SRC.items():
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
    Ds, Gs = [], []
    for fn in sorted(os.listdir(P)):
        if not fn.endswith('.jsonl.gz'): continue
        for line in gzip.open(os.path.join(P, fn), 'rt'):
            r = json.loads(line); Ds.append(r['D']); Gs.append(r['g'])
    for nm, s, a in (('kl', 1, 0), ('zap', -1, 1)):
        idx = [i for i, d in enumerate(Ds) if (d > 0) == (s > 0) and d in lv]
        Dv = np.array([Ds[i] for i in idx], np.int64); L = np.array([lv[Ds[i]] for i in idx]); C = np.zeros((len(idx), nb))
        for j, i in enumerate(idx):
            t = np.array(Gs[i], float)
            if len(t) == 0: continue
            x = ((t / 2) * math.log(abs(Ds[i]) / math.pi) + np.imag(lngamma_c((0.5 + a + 1j * t) / 2))) / math.pi
            C[j] = np.histogram(x, edges)[0]
        def stat(wt): return ((wt * L) @ C / np.sum(wt * L) - wt @ C / np.sum(wt)) / bw, (wt * L) @ C / np.sum(wt * L) / bw, wt @ C / np.sum(wt) / bw
        d, dw, d0 = stat(np.ones(len(L))); B = np.array([stat(rng.poisson(1.0, len(L)).astype(float))[0] for _ in range(300)]); e = B.std(axis=0, ddof=1)
        key = '%s_%s' % (F, nm); Nx = H2['skupiny'][key]['N_x1']; du, (u0, u1) = usp_N(Nx)
        chi2 = float(np.sum(((d - du) / e) ** 2)); p = float(stats.chi2.sf(chi2, nb))
        cov = np.cov(B.T); r = d - du; chi2c = float(r @ np.linalg.solve(cov, r))
        lim = 1 + np.sin(2*np.pi*(edges[:-1]+edges[1:])/2)/(2*np.pi*(edges[:-1]+edges[1:])/2) - 2*np.sin(np.pi*(edges[:-1]+edges[1:])/2)**2/(np.pi*(edges[:-1]+edges[1:])/2)**2 - (1 - np.sin(2*np.pi*(edges[:-1]+edges[1:])/2)/(2*np.pi*(edges[:-1]+edges[1:])/2))
        sh = stat(np.ones(len(L)))  # nulový model: promíchání L
        Lp = rng.permutation(L); dnul = ((Lp @ C) / Lp.sum() - C.sum(0) / len(L)) / bw
        out['skupiny'][key] = dict(n=len(L), N_x1=Nx, rozdil=d.tolist(), chyba=e.tolist(), usp=du.tolist(), vazena=dw.tolist(), nevazena=d0.tolist(), usp_vazena=u1.tolist(), usp_nevazena=u0.tolist(), limita=lim.tolist(), chi2=chi2, p=p, chi2_kov=chi2c, nul_max=float(np.max(np.abs(dnul / e))))
        print('%-9s n %6d N_x1 %.2f  χ² %.1f/15 p %.2e (s kovariancí %.1f)  nulový max|z| %.1f  rozdíl koš 0–0,2: data %+.3f USp %+.3f limita %+.3f; koš 0,4–0,6: %+.3f / %+.3f' % (key, len(L), Nx, chi2, p, chi2c, out['skupiny'][key]['nul_max'], d[0], du[0], lim[0], d[2], du[2]), flush=True)
ps = [g['p'] for g in out['skupiny'].values()]
out['verdikt'] = 'O' if all(p > 1e-3 for p in ps) else ('Z' if sum(p < 1e-4 for p in ps) >= 2 else 'N'); print('verdikt', out['verdikt'])
json.dump(out, open('/home/claude/br78/data/h3_hust.json', 'w'), indent=1)
