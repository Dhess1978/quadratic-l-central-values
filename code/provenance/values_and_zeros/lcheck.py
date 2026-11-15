"""Independent check of L(1/2, chi_D): smoothed approximate functional equation,
L(1/2) = 2 sum_n chi_D(n) n^(-1/2) Gamma(c, pi n^2/|D|)/Gamma(c), c = 1/4 (D>0) or 3/4 (D<0).
chi_D(n) by own Kronecker symbol (completely multiplicative fill over primes). No PARI."""
import numpy as np, math
from scipy.special import gammaincc
def spf_table(N):
    spf = np.zeros(N + 1, dtype=np.int64)
    for i in range(2, N + 1):
        if spf[i] == 0:
            spf[i::i][spf[i::i] == 0] = i
    return spf
def kron_p(D, p):
    if p == 2:
        if D % 2 == 0: return 0
        return 1 if D % 8 in (1, 7) else -1
    r = D % p
    if r == 0: return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1
def chi_vec(D, N, spf):
    chi = np.zeros(N + 1, dtype=np.int64); chi[1] = 1
    cache = {}
    for n in range(2, N + 1):
        p = int(spf[n])
        if p not in cache: cache[p] = kron_p(D, p)
        chi[n] = cache[p] * chi[n // p]
    return chi
def L_half(D, spf, cut=46.0):
    q = abs(D); c = 0.25 if D > 0 else 0.75
    N = int(math.sqrt(cut * q / math.pi)) + 2
    chi = chi_vec(D, N, spf)
    n = np.arange(1, N + 1, dtype=float)
    w = gammaincc(c, math.pi * n * n / q)   # regularized upper incomplete gamma
    return 2.0 * float(np.sum(chi[1:] * w / np.sqrt(n)))
def L_half_mp(D, dps=30, cut=60):
    import mpmath as mp
    mp.mp.dps = dps
    q = abs(D); c = mp.mpf(1)/4 if D > 0 else mp.mpf(3)/4
    N = int(math.sqrt(cut * q / math.pi)) + 2
    s = mp.mpf(0); g = mp.gamma(c)
    for n in range(1, N + 1):
        x = 1
        m = n; ch = 1
        # Kronecker via factorization by trial (small n)
        k = 2
        while m > 1 and ch != 0:
            if k * k > m: ch *= kron_p(D, m); m = 1; break
            while m % k == 0: ch *= kron_p(D, k); m //= k
            k += 1
        if ch == 0: continue
        s += ch * mp.gammainc(c, mp.pi * n * n / q) / (g * mp.sqrt(n))
    return 2 * s
