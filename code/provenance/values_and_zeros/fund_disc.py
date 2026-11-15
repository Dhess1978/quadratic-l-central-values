"""All fundamental discriminants with 3 <= |D| <= X (both signs), via a squarefree sieve."""
import numpy as np
def squarefree_mask(N):
    m = np.ones(N + 1, dtype=bool); m[0] = False
    p = 2
    lim = int(N ** 0.5) + 1
    sieve = np.ones(lim + 1, dtype=bool); sieve[:2] = False
    for i in range(2, int(lim ** 0.5) + 1):
        if sieve[i]: sieve[i*i::i] = False
    for p in np.nonzero(sieve)[0]:
        m[p*p::p*p] = False
    return m
def fundamental(X):
    sf = squarefree_mask(X)
    n = np.arange(X + 1)
    out = []
    for s in (1, -1):
        D = s * n  # candidate D = s*n
        # odd case: D squarefree, D = 1 mod 4
        odd = sf & ((D % 4) == 1) & (n >= 2)
        out.append(D[odd])
        # even case: D = 4m, m squarefree, m = 2 or 3 mod 4 (m = D/4 with sign)
        m = n // 4
        ok = (n % 4 == 0) & (n >= 4)
        mm = s * m
        cond = ok & sf[m] & (((mm % 4) == 2) | ((mm % 4) == 3))
        out.append(D[cond])
    return np.sort(np.concatenate(out))
if __name__ == '__main__':
    import sys
    X = int(float(sys.argv[1]))
    F = fundamental(X)
    print(X, len(F), (F > 0).sum(), (F < 0).sum(), F[:6], F[-3:])
