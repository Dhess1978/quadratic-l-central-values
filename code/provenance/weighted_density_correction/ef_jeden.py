"""ověření konvencí explicitního vzorce na jednotlivých L(s, χ_D): nuly z PARI do u = 40 proti prvočíselné straně"""
import math, numpy as np, cypari2, scipy.special as sp
from scipy import integrate
pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9)
phi = lambda u: np.sinc(u / 2) ** 4
def M4(x):
    x = np.abs(x); return np.where(x <= 1, 2 / 3 - x ** 2 + x ** 3 / 2, np.where(x <= 2, (2 - x) ** 3 / 6, 0.0))
phihat = lambda nu: 2 * M4(2 * nu)
for nu in (0, 0.3, 0.7, 0.95):
    v = integrate.quad(lambda u: phi(u) * math.cos(2 * math.pi * u * nu), -200, 200, limit=4000)[0]; print('φ̂(%.2f) numericky %.6f vzorec %.6f' % (nu, v, phihat(nu)))
def arch(l, a):
    f = lambda u: phi(u) * np.real(sp.digamma((0.5 + a + 2j * math.pi * u / l) / 2))
    return 4 / 3 + integrate.quad(f, -400, 400, limit=8000, points=[0])[0] / l
def kron(D, n): return int(pari.kronecker(D, n))
for D in (1000001, -1000003, 10001, -10003):
    a = 0 if D > 0 else 1; l = math.log(abs(D) / math.pi); T = 40 * 2 * math.pi / l
    g = np.array([float(v) for v in pari.lfunzeros(pari.lfuncreate(D), T, 32) if float(v) > 1e-9])
    Z = 2 * phi(g * l / (2 * math.pi)).sum()
    S = 0.0; N = int(abs(D) / math.pi)
    for p in pari.primes([2, N]):
        p = int(p); pk = p; k = 1
        while pk <= N:
            c = kron(D, p) ** k
            if c: S += math.log(p) * c / math.sqrt(pk) * phihat(math.log(pk) / l)
            pk *= p; k += 1
    P = arch(l, a) - 2 / l * S
    print('D %8d  nul %4d  strana nul %.8f  strana prvočísel %.8f  rozdíl %.1e' % (D, len(g), Z, P, Z - P), flush=True)
