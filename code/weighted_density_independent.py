"""Study V — independent implementation of the prediction of the weighted (L(1/2)) and unweighted one-level density of low zeros of quadratic
Dirichlet L-functions from the ratios conjecture (Conrey-Farmer-Zirnbauer recipe, Fazzari's weighted density), written without the code of the
first implementation (code/provenance/weighted_density/fazzari.py): local factors as numerically summed series (complete homogeneous polynomials
h_e(x, y) by the recursion h_e = x h_{e-1} + y^e), Euler product up to 5e5, derivative in beta by a central difference (1e-5), limit alpha -> 0 by
Richardson extrapolation from alpha = +-h, +-2h (h = 1e-3), family average <chi_D(k)> = prod_{p | k} p/(p+1) for squares k.
It found the error of the first implementation (Section 4.6 of the article) and reproduces the corrected prediction published in
data/results/study_v_results.json (weighted_density.corrected.groups) for every family whose discriminants are in data/central_values
(B10, C100, S7, D100, D100+D100b, D100c, all three at 1e8, S9L at 1e9, E; R10 is a random subset and is skipped).
usage: python3 -I weighted_density_independent.py REPOSITORY_ROOT [workers]      (about 15 minutes with 2 workers)"""
import sys, time, math, numpy as np, mpmath as mp
from multiprocessing import Pool
import scipy.special as sp
mp.mp.dps = 25
PM = 500000; s = np.ones(PM + 1, bool); s[:2] = False
for q in range(2, int(PM ** 0.5) + 1):
    if s[q]: s[q * q::q] = False
P = np.nonzero(s)[0].astype(float); LP = np.log(P); R = P / (P + 1)
MAL = P < 60; VEL = ~MAL; EMAX_MAL, EMAX_VEL = 120, 26
def mocn(e): return np.exp(-e * LP)     # p^{-e}, complex e
def _rady(x, y, emax):
    """sums over even e >= 2 and over odd e of h_e(x, y), numerically"""
    h = np.ones_like(x); ye = np.ones_like(y); Ss = np.zeros_like(x); So = np.zeros_like(x)
    for e in range(1, emax + 1):
        ye = ye * y; h = x * h + ye
        if e % 2 == 0: Ss = Ss + h
        else: So = So + h
    return Ss, So
def logA3(al, be, ga):
    x, y, z = mocn(0.5 + al), mocn(0.5 + be), mocn(0.5 + ga); tot = 0j
    for m, em in ((MAL, EMAX_MAL), (VEL, EMAX_VEL)):
        xm, ym, zm, rm = x[m], y[m], z[m], R[m]; Ss, So = _rady(xm, ym, em)
        G = 1 + rm * Ss - zm * rm * So
        tot += np.sum(np.log(G) + np.log(1 - xm * xm) + np.log(1 - ym * ym) + np.log(1 - xm * ym) - np.log(1 - xm * zm) - np.log(1 - ym * zm))
    return tot
def logA2(be, ga):
    y, z = mocn(0.5 + be), mocn(0.5 + ga); tot = 0j
    for m, em in ((MAL, EMAX_MAL), (VEL, EMAX_VEL)):
        ym, zm, rm = y[m], z[m], R[m]; Ss, So = _rady(np.zeros_like(ym), ym, em)   # x = 0: h_e = y^e
        G = 1 + rm * Ss - zm * rm * So
        tot += np.sum(np.log(G) + np.log(1 - ym * ym) - np.log(1 - ym * zm))
    return tot
def logA1(al):
    x = mocn(0.5 + al); tot = 0j
    for m, em in ((MAL, EMAX_MAL), (VEL, EMAX_VEL)):
        xm, rm = x[m], R[m]; Ss, _ = _rady(np.zeros_like(xm), xm, em)
        tot += np.sum(np.log(1 + rm * Ss) + np.log(1 - xm * xm))
    return tot
def ze(z): return complex(mp.zeta(mp.mpc(z.real, z.imag)))
def F3(al, be, ga):
    Y = ze(1 + 2 * al) * ze(1 + 2 * be) * ze(1 + al + be) / (ze(1 + al + ga) * ze(1 + be + ga)); return Y * np.exp(logA3(al, be, ga))
def F2(be, ga): return ze(1 + 2 * be) / ze(1 + be + ga) * np.exp(logA2(be, ga))
DL = 1e-5; H = 1e-3; AL = [H, -H, 2 * H, -2 * H]
def kus_t(t):
    g = 1j * t; out = np.zeros((4, 4), complex)    # [alpha index, swap]
    for i, al in enumerate(AL):
        for k, (sa, sb) in enumerate(((1, 1), (-1, 1), (1, -1), (-1, -1))):
            f = lambda b: F3(sa * al, sb * b, g)
            out[i, k] = (f(g + DL) - f(g - DL)) / (2 * DL)
    u = np.array([(F2(g + DL, g) - F2(g - DL, g)) / (2 * DL), (F2(-(g + DL), g) - F2(-(g - DL), g)) / (2 * DL)])
    return t, out, u
def X(s, l, a):   # X_D(s) = (|D|/π)^{1/2 - s} Γ((1 - s + a)/2) / Γ((s + a)/2)
    return np.exp(-(s - 0.5) * l + sp.loggamma((1 - s + a) / 2) - sp.loggamma((s + a) / 2))
def rich(g): return (2 * (g[..., 0] + g[..., 1]) - 0.5 * (g[..., 2] + g[..., 3])) / 3     # α = h, −h, 2h, −2h → α = 0
def hustoty(l, a):
    global TG, W, U, Mk, DT
    Xa = X(0.5 + np.array(AL), l, a); Xt = X(0.5 + 1j * TG, l, a)
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
    import os, math, json
    ROOT = sys.argv[1]; NW = int(sys.argv[2]) if len(sys.argv) > 2 else 2; t0 = time.time()
    DT = 0.0025; TG = (np.arange(int(3.5 / DT)) + 0.5) * DT
    Mk = np.array([[ze(1 + 2 * al) * np.exp(logA1(al)), ze(1 - 2 * al) * np.exp(logA1(-al))] for al in AL])
    W = np.zeros((len(TG), 4, 4), complex); U = np.zeros((len(TG), 2), complex)
    with Pool(NW) as pool:
        for n, (t, o, u) in enumerate(pool.imap(kus_t, TG, chunksize=8)):
            j = int(round(t / DT - 0.5)); W[j] = o; U[j] = u
            if n % 200 == 0: print('pieces', n, len(TG), round(time.time() - t0), 's', flush=True)
    R = json.load(open(os.path.join(ROOT, 'data/results/study_v_results.json')))['weighted_density']['corrected']['groups']
    CV = os.path.join(ROOT, 'data/central_values'); ld = lambda n: np.load(os.path.join(CV, n + '.npz'))['D']
    V = np.concatenate([ld(f[:-4]) for f in sorted(os.listdir(CV)) if f.startswith('V_part')])
    S = {'B10': V[(np.abs(V) > 10 ** 4) & (np.abs(V) <= 10 ** 5)], 'C100': V[(np.abs(V) > 10 ** 6) & (np.abs(V) <= 2 * 10 ** 6)], 'S7': V[(np.abs(V) > 7500000) & (np.abs(V) <= 10 ** 7)], 'D100': ld('D100'), 'D100_D100b': np.concatenate([ld('D100'), ld('D100b')]),
         'D100c': ld('D100c'), 'D100_D100b_D100c': np.concatenate([ld('D100'), ld('D100b'), ld('D100c')]), 'S9L': ld('S9L'), 'E': ld('E10_E100')}
    mx = 0.0
    for F, Ds in S.items():
        for zn, sg, a in (('even', 1, 0), ('odd', -1, 1)):
            g = R[F + '_' + zn]; lv = np.log(np.abs(Ds[np.sign(Ds) == sg]) / math.pi); pr, vz, nv = predikce(lv, a)
            d = np.max(np.abs(pr - np.array(g['prediction'])) / np.array(g['error'])); mx = max(mx, d)
            c2 = float(np.sum(((np.array(g['data']) - pr) / np.array(g['error'])) ** 2))
            print('%-18s %-4s n %7d  |independent - published| <= %.4f errors   chi2 %.2f (published %.2f)' % (F, zn, len(lv), d, c2, g['chi2']), flush=True)
    print('largest difference %.4f errors (%s); %d s' % (mx, 'PASS' if mx < 0.05 else 'FAIL', time.time() - t0))
