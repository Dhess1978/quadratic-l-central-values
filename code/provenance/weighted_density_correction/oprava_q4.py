"""ověření příčiny: v fazzari.py (34.79) má čtvrtý člen (dvojí výměna α → −α, β → −β) Γ-podíl gS(½+α+β) místo gS(½+α)·gS(½+β).
Opraví se jen ten faktor v uložených kusech a přepočte předpověď stejným kódem jako 34.79 (vahova, nevazena, kose)."""
import json, math, numpy as np, mpmath as mp, scipy.special as sp
K = np.load('/home/claude/br79/data/fazzari_kusy.npz'); TG = K['TG']; AJ = K['AJ']; DT = 0.002; NT = len(TG)
def gS(z, a): return complex(mp.gamma(mp.mpc((1 - z.real + a) / 2, -z.imag / 2)) / mp.gamma(mp.mpc((z.real + a) / 2, z.imag / 2)))
PC = {}
for a in (0, 1):
    Q = K['Q%d' % a].copy()
    for i, t in enumerate(TG):
        D = 1j * float(t)
        for j, A in enumerate(AJ): Q[i, j, 3] *= gS(0.5 + A, a) * gS(0.5 + D, a) / gS(0.5 + A + D, a)
    PC[a] = (Q, K['M%d' % a])
IMLG = {k: np.imag(sp.loggamma(0.25 + 0.5 * k + 0.5j * TG)) for k in (0, 1)}
def vahova(l, a, PC):
    Q, Mk = PC[a]; M = np.mean(Mk[:, 0] + np.exp(-AJ * l) * Mk[:, 1])
    E = np.stack([np.ones((NT, len(AJ))), np.exp(-AJ * l)[None, :].repeat(NT, 0), np.exp(-1j * TG * l)[:, None].repeat(len(AJ), 1), np.exp(-(AJ[None, :] + 1j * TG[:, None]) * l)], -1)
    S = np.mean(np.sum(Q * E, -1), 1); return M, np.real(S) / math.pi
def kose(l, a, c, edges):
    tt = np.concatenate([[0], TG]); c0 = 3 * c[0] - 3 * c[1] + c[2]; cc = np.concatenate([[c0], c]); C = np.concatenate([[0], np.cumsum(0.5 * (cc[1:] + cc[:-1]) * DT)])
    th = (0.5 * tt * l + np.concatenate([[0], IMLG[a]])) / math.pi; te = np.interp(edges, th, tt); return np.diff(np.interp(te, tt, C))
EDG = np.linspace(0, 1, 6); BW = 0.2
PC0 = {a: (K['Q%d' % a], K['M%d' % a]) for a in (0, 1)}
for lc in (9.62, 13.06, 17.29, 21.88):
    for a in (0, 1):
        M1, c1 = vahova(lc, a, PC0); M2, c2 = vahova(lc, a, PC)
        v1 = (M1.real * BW + kose(lc, a, c1, EDG)) / (M1.real * BW); v2 = (M2.real * BW + kose(lc, a, c2, EDG)) / (M2.real * BW)
        print('l %.2f a %d  vážená 34.79 %s  opravená %s' % (lc, a, np.round(v1, 4), np.round(v2, 4)))
