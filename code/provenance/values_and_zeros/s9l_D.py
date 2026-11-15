"""34.85b: seznam základních diskriminantů rodiny S9L (10⁹ < |D| ≤ 10⁹ + 8·10⁵, obě znaménka) vlastním sítem v okně. usage: python3 -I s9l_D.py VYSTUP.npz"""
import sys, math, numpy as np
LO, HI = 10 ** 9, 10 ** 9 + 800000
def bezctvercova(a, b):      # pole bool pro n = a … b
    n = np.arange(a, b + 1); ok = np.ones(len(n), bool)
    for p in range(2, int(math.isqrt(b)) + 1):
        if p > 3 and (p % 2 == 0 or p % 3 == 0): continue
        q = p * p; s = (-a) % q; ok[s::q] = False
    return n, ok
n, sf = bezctvercova(LO + 1, HI); kl = [n[sf & (n % 4 == 1)]]; zap = [-n[sf & (n % 4 == 3)]]
m, sm = bezctvercova(LO // 4 + 1, HI // 4); ok = sm & (4 * m > LO); kl.append(4 * m[ok & ((m % 4 == 2) | (m % 4 == 3))]); zap.append(-4 * m[ok & ((m % 4 == 1) | (m % 4 == 2))])
D = np.concatenate(kl + zap).astype(np.int64); D = D[np.argsort(np.abs(D), kind='stable')]
np.savez(sys.argv[1], D=D); print('S9L základních D', len(D), 'sudých', int((D > 0).sum()), 'lichých', int((D < 0).sum()), 'min', int(np.abs(D).min()), 'max', int(np.abs(D).max()))
