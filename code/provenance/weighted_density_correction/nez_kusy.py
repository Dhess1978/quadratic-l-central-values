"""34.85 odv-lv-nez: NEZÁVISLÁ implementace předpovědi domněnky poměrů pro váženou (L(½)) a neváženou hustotu nízkých nul kvadratické rodiny.
Nesdílí kód s fazzari.py: lokální faktory jako numericky sčítané řady (úplné homogenní polynomy h_e(x, y), rekurze h_e = x h_{e-1} + y^e),
Eulerův součin do 5·10^5, derivace podle β centrální diferencí (δ = 1e-5), limita α → 0 Richardsonem z α = ±h, ±2h (h = 1e-3),
recept CFZ: ⟨L(½+α)L(½+β)/L(½+γ)⟩ = Σ přes výměny α → −α, β → −β (s X(½+α), X(½+β)) z Y·A, ⟨χ_D(k)⟩ = Π_{p|k} p/(p+1) pro čtverce k.
Výstup: kusy nezávislé na D (pro sudé i liché jsou stejné; znaménko a vodič vstupují jen přes X) na mřížce t."""
import sys, time, math, numpy as np, mpmath as mp
from multiprocessing import Pool
mp.mp.dps = 25
PM = 500000; s = np.ones(PM + 1, bool); s[:2] = False
for q in range(2, int(PM ** 0.5) + 1):
    if s[q]: s[q * q::q] = False
P = np.nonzero(s)[0].astype(float); LP = np.log(P); R = P / (P + 1)
MAL = P < 60; VEL = ~MAL; EMAX_MAL, EMAX_VEL = 120, 26
def mocn(e): return np.exp(-e * LP)     # p^{-e}, e komplexní
def _rady(x, y, emax):
    """Σ_{e sudé ≥ 2} h_e(x, y) a Σ_{e liché} h_e(x, y) numericky"""
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
    g = 1j * t; out = np.zeros((4, 4), complex)    # [α index, výměna]
    for i, al in enumerate(AL):
        for k, (sa, sb) in enumerate(((1, 1), (-1, 1), (1, -1), (-1, -1))):
            f = lambda b: F3(sa * al, sb * b, g)
            out[i, k] = (f(g + DL) - f(g - DL)) / (2 * DL)
    u = np.array([(F2(g + DL, g) - F2(g - DL, g)) / (2 * DL), (F2(-(g + DL), g) - F2(-(g - DL), g)) / (2 * DL)])
    return t, out, u
if __name__ == '__main__':
    t0 = time.time(); DT = 0.0025; TG = (np.arange(int(3.5 / DT)) + 0.5) * DT
    Mk = np.array([[ze(1 + 2 * al) * np.exp(logA1(al)), ze(1 - 2 * al) * np.exp(logA1(-al))] for al in AL])
    W = np.zeros((len(TG), 4, 4), complex); U = np.zeros((len(TG), 2), complex)
    with Pool(2) as pool:
        for n, (t, o, u) in enumerate(pool.imap(kus_t, TG, chunksize=8)):
            j = int(round(t / DT - 0.5)); W[j] = o; U[j] = u
            if n % 100 == 0: print(n, len(TG), round(time.time() - t0), 's', flush=True)
    np.savez('/home/claude/br85/data/nez_kusy.npz', TG=TG, W=W, U=U, Mk=Mk, AL=np.array(AL), H=H, DL=DL, PM=PM)
    print('hotovo', round(time.time() - t0), 's', flush=True)
