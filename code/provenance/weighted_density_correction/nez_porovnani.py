"""34.85 odv-lv-nez: předpověď z nezávislé implementace (nez_kusy.npz) pro koše 0,2 v x ∈ [0, 1] a porovnání s dosavadní (fazzari.py / lu*.py) v jednotkách chyby dat."""
import json, math, numpy as np, scipy.special as sp
K = np.load('/home/claude/br85/data/nez_kusy.npz'); TG, W, U, Mk, AL, DT = K['TG'], K['W'], K['U'], K['Mk'], K['AL'], float(K['TG'][1] - K['TG'][0])
def X(s, l, a):   # X_D(s) = (|D|/π)^{1/2 - s} Γ((1 - s + a)/2) / Γ((s + a)/2)
    return np.exp(-(s - 0.5) * l + sp.loggamma((1 - s + a) / 2) - sp.loggamma((s + a) / 2))
def rich(g): return (2 * (g[..., 0] + g[..., 1]) - 0.5 * (g[..., 2] + g[..., 3])) / 3     # α = h, −h, 2h, −2h → α = 0
def hustoty(l, a):
    Xa = X(0.5 + AL, l, a); Xt = X(0.5 + 1j * TG, l, a)
    M = rich(Mk[:, 0] + Xa * Mk[:, 1])
    S = rich(W[:, :, 0] + Xa[None, :] * W[:, :, 1] + Xt[:, None] * W[:, :, 2] + (Xa[None, :] * Xt[:, None]) * W[:, :, 3])
    Su = U[:, 0] + Xt * U[:, 1]
    th = 0.5 * TG * l + np.imag(sp.loggamma((0.5 + a + 1j * TG) / 2)); thp = 0.5 * l + 0.5 * np.real(sp.digamma((0.5 + a + 1j * TG) / 2))
    return M.real, thp / math.pi, np.real(S) / math.pi, np.real(Su) / math.pi
def kose(l, a, edges=np.linspace(0, 1, 6)):
    M, rho0, s_w, s_u = hustoty(l, a)
    te = np.concatenate([[0], (np.arange(len(TG)) + 1) * DT]); th = 0.5 * te * l + np.imag(sp.loggamma((0.5 + a + 1j * te) / 2)); x = th / math.pi
    tk = np.interp(edges, x, te)
    def integ(f): C = np.concatenate([[0], np.cumsum(f) * DT]); return np.diff(np.interp(tk, te, C))
    return M, integ(rho0), integ(s_w), integ(s_u)
def predikce(lv, a, nb=40):
    ed = np.linspace(lv.min(), lv.max() + 1e-9, nb + 1); idx = np.digitize(lv, ed) - 1
    num = np.zeros(5); den = 0.0; nu = np.zeros(5); nw = 0.0; BW = 0.2
    for i in range(nb):
        m = idx == i
        if not m.any(): continue
        w = m.sum(); lc = lv[m].mean(); M, r0, sw, su = kose(lc, a)
        num += w * (M * r0 + sw); den += w * M * BW; nu += w * (r0 + su); nw += w * BW
    return num / den - nu / nw, num / den, nu / nw
if __name__ == '__main__':
    R = json.load(open('/home/claude/pub5/repo/data/results/study_v_results.json'))['weighted_density']
    CV = '/home/claude/pub5/repo/data/central_values/'
    def ld(n): z = np.load(CV + n + '.npz'); return z['D']
    V = np.concatenate([np.load(CV + f)['D'] for f in sorted(__import__('os').listdir(CV)) if f.startswith('V_part')])
    Dsets = {'B10': V[(np.abs(V) > 10 ** 4) & (np.abs(V) <= 10 ** 5)], 'C100': V[(np.abs(V) > 10 ** 6) & (np.abs(V) <= 2 * 10 ** 6)], 'R10': np.load('/home/claude/br76/data/R10.npz')['D'],
             'D100': ld('D100'), 'D100_D100b': np.concatenate([ld('D100'), ld('D100b')]), 'D100c': ld('D100c'), 'D100_D100b_D100c': np.concatenate([ld('D100'), ld('D100b'), ld('D100c')]), 'E': ld('E10_E100')}
    out = {}; mx = 0.0
    for F, Ds in Dsets.items():
        for zn, sg, a in (('even', 1, 0), ('odd', -1, 1)):
            g = R['groups'][F + '_' + zn] if F in ('B10', 'C100', 'R10', 'D100') else R[F][zn]
            Dz = Ds[np.sign(Ds) == sg]; lv = np.log(np.abs(Dz) / math.pi)
            pr, vz, nv = predikce(lv, a); old = np.array(g['pred']); e = np.array(g['chyba']); d = np.array(g['rozdil'])
            r = np.abs(pr - old) / e; mx = max(mx, r.max())
            c2n = float(np.sum(((d - pr) / e) ** 2)); c2o = float(np.sum(((d - old) / e) ** 2))
            out[F + '_' + zn] = dict(n=int(len(Dz)), n_old=g.get('n'), pred_nez=pr.tolist(), pred_old=old.tolist(), max_rozdil_v_chybach=float(r.max()), chi2_nez=c2n, chi2_old=c2o, vazena_nez=vz.tolist(), nevazena_nez=nv.tolist())
            print('%-22s n %7d (%s)  max|nová − dosavadní|/chyba %.4f   χ² nová %.1f  dosavadní %.1f' % (F + '_' + zn, len(Dz), g.get('n'), r.max(), c2n, c2o), flush=True)
    out['max'] = mx; out['verdikt'] = 'ODVOZENO' if mx <= 0.05 else ('ZAMÍTNUTO' if mx > 0.5 else 'NELZE ODVODIT')
    print('největší rozdíl %.4f chyby → %s' % (mx, out['verdikt']))
    json.dump(out, open('/home/claude/br85/data/nez_porovnani.json', 'w'), indent=1)
