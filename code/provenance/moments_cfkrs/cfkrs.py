"""CFKRS moment polynomials Q_k(x), x = ln|d|, for L(1/2, chi_d) over fundamental discriminants of one sign.
Q_k(x) = ((-1)^{k(k-1)/2} 2^k / k!) (2 pi i)^{-k} oint G(z) Delta(z^2)^2 / prod z_j^{2k-1} exp(x/2 sum z) dz,
G(z) = A_k(z) prod_j X(1/2+z_j, a)^{-1/2} prod_{i<=j} zeta(1+z_i+z_j),
X(1/2+z, a) = pi^z Gamma((1/2-z+a)/2)/Gamma((1/2+z+a)/2), a = 0 (d>0), 1 (d<0),
A_k(z) = prod_p prod_{i<=j}(1-p^{-1-z_i-z_j}) * (1/(p+1) + p/(p+1) * (prod(1-x_j)^-1 + prod(1+x_j)^-1)/2), x_j = p^{-1/2-z_j}.
Numerics: trapezoid on circles (different radius per variable), primes <= P exact, tail by integral approximation."""
import numpy as np, itertools, math
from scipy.special import loggamma
import mpmath as mp
_STI = [float(mp.stieltjes(n)) for n in range(45)]
def zeta1p(w):  # zeta(1+w) for complex array w, |w| <= ~1.2
    s = 1.0 / w
    term = np.ones_like(w)
    for n in range(45):
        s = s + ((-1) ** n) * _STI[n] * term / math.factorial(n)
        term = term * w
    return s
def primes_upto(P):
    s = np.ones(P + 1, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def _deg4(k):
    """Coefficients of the degree-4 part of log A_p in x_j = p^{-1/2-z_j}: h4 - h2^2/2 - (1/2) sum_{i<=j} x_i^2 x_j^2 (34.76 tail fix)."""
    from collections import defaultdict
    c = defaultdict(float); pairs = [(i, j) for i in range(k) for j in range(i, k)]
    for a in itertools.combinations_with_replacement(range(k), 4):
        al = tuple(a.count(j) for j in range(k)); c[al] += 1.0                      # h4
    for (i, j) in pairs:
        for (l, m) in pairs:
            al = [0] * k; al[i] += 1; al[j] += 1; al[l] += 1; al[m] += 1; c[tuple(al)] -= 0.5   # -h2^2/2
        al = [0] * k; al[i] += 2; al[j] += 2; c[tuple(al)] -= 0.5                   # -(1/2) sum x_i^2 x_j^2
    return [(np.array(a, float), v) for a, v in c.items() if v != 0]
def logA(Z, P=100000, tail='v2'):
    """Z: array (npts, k) complex. returns log A_k(Z). tail='v2' (34.76): all p^-2 terms of log A_p for p > P
    (degree-4 monomials and -h2/(p+1)); tail='v1' (10. 11. dopoledne): only -sum_{i<=j} p^{-2-z_i-z_j} (at z = 0: -k(k+1)/2 instead of e.g. -30 for k = 4)."""
    npts, k = Z.shape
    out = np.zeros(npts, dtype=complex)
    pairs = [(i, j) for i in range(k) for j in range(i, k)]
    for p in primes_upto(P):
        lp = math.log(p)
        x = np.exp(-(0.5 + Z) * lp)          # p^{-1/2-z_j}
        prod_m = np.prod(1.0 / (1.0 - x), axis=1); prod_p = np.prod(1.0 / (1.0 + x), axis=1)
        loc = 1.0 / (p + 1) + p / (p + 1.0) * 0.5 * (prod_m + prod_p)
        for i, j in pairs:
            loc = loc * (1.0 - x[:, i] * x[:, j])
        out += np.log(loc)
    # tail p > P: local log ~ -sum_{i<=j} p^{-2-z_i-z_j} * (1 + ...) ; leading term only
    # sum_{p>P} p^{-2-w} ~ P^{-1-w} / ((1+w) ln P)
    lP = math.log(P)
    from scipy.special import exp1
    T = (lambda w: exp1((1 + w) * lP)) if tail == 'v2' else (lambda w: np.exp(-(1 + w) * lP) / ((1 + w) * lP))   # v2: sum_{p>P} p^{-2-w} ~ int dt/(t^{2+w} ln t) = E1((1+w) ln P)
    for i, j in pairs:
        out += -T(Z[:, i] + Z[:, j])
    if tail == 'v2':
        for al, c in _deg4(k):
            out += c * T(Z @ al)
    return out
def Qk_coeffs(k, a, N=24, radii=None, P=100000, tail='v2'):
    """Return polynomial coefficients c_m (m = 0..k(k+1)/2) with Q_k(x) = sum c_m x^m."""
    if radii is None: radii = [0.22 + 0.031 * j for j in range(k)]
    th = 2 * np.pi * (np.arange(N) + 0.5) / N
    grids = [r * np.exp(1j * (th + 0.37 * j)) for j, r in enumerate(radii)]
    Z = np.array(list(itertools.product(*grids)))            # (N^k, k)
    # dz_j/(2 pi i) on circle = z_j dtheta/(2 pi) -> weight z_j / N
    w = np.prod(Z, axis=1) / N ** k
    # G(z)
    lg = logA(Z, P, tail)
    for j in range(k):
        z = Z[:, j]; c = (0.5 + a) / 2
        lX = z * math.log(math.pi) + loggamma(c - z / 2) - loggamma(c + z / 2)
        lg = lg - 0.5 * lX
    G = np.exp(lg)
    for i in range(k):
        for j in range(i, k):
            G = G * zeta1p(Z[:, i] + Z[:, j])
    V = np.ones(len(Z), dtype=complex)
    for i in range(k):
        for j in range(i + 1, k):
            V = V * (Z[:, i] ** 2 - Z[:, j] ** 2) ** 2
    pre = ((-1) ** (k * (k - 1) // 2)) * 2 ** k / math.factorial(k)
    R = pre * G * V / np.prod(Z ** (2 * k - 1), axis=1) * w
    s = Z.sum(axis=1) / 2
    deg = k * (k + 1) // 2
    coeffs = []
    sm = np.ones_like(s)
    for m in range(deg + 1):
        coeffs.append(np.sum(R * sm).real / math.factorial(m))
        sm = sm * s
    return np.array(coeffs)
def a_k(k, P=100000, tail='v2'):
    return float(np.exp(logA(np.zeros((1, k), dtype=complex), P, tail)[0].real))
def ks_leading(k):
    c = 1.0
    for j in range(1, k + 1): c *= math.factorial(j) / math.factorial(2 * j)
    return c
