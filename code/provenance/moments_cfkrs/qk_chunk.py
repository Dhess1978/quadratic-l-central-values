"""34.77: CFKRS Q_k po částech mřížky (první dvě souřadnice, N² částí) s kontrolním bodem — po pádu pokračuje.
Použití: python3 qk_chunk.py k a N r1,r2,... VYSTUP.json ADRESAR_CHECKPOINTU"""
import sys, os, json, time, math, itertools, numpy as np
sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
from cfkrs import logA, zeta1p
from scipy.special import loggamma
k, a, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); radii = [float(x) for x in sys.argv[4].split(',')]; OUT, CK = sys.argv[5], sys.argv[6]
os.makedirs(CK, exist_ok=True); P = 100000; tail = 'v2'; deg = k * (k + 1) // 2
th = 2 * np.pi * (np.arange(N) + 0.5) / N
grids = [r * np.exp(1j * (th + 0.37 * j)) for j, r in enumerate(radii)]
t0 = time.time()
for i0 in range(N * N):
    f = os.path.join(CK, 'c%03d.json' % i0)
    if os.path.exists(f): continue
    t = time.time()
    Z = np.array([(grids[0][i0 // N], grids[1][i0 % N]) + p for p in itertools.product(*grids[2:])])
    w = np.prod(Z, axis=1) / N ** k
    lg = logA(Z, P, tail)
    for j in range(k):
        z = Z[:, j]; c = (0.5 + a) / 2
        lg = lg - 0.5 * (z * math.log(math.pi) + loggamma(c - z / 2) - loggamma(c + z / 2))
    G = np.exp(lg)
    for i in range(k):
        for j in range(i, k): G = G * zeta1p(Z[:, i] + Z[:, j])
    V = np.ones(len(Z), dtype=complex)
    for i in range(k):
        for j in range(i + 1, k): V = V * (Z[:, i] ** 2 - Z[:, j] ** 2) ** 2
    R = ((-1) ** (k * (k - 1) // 2)) * 2 ** k / math.factorial(k) * G * V / np.prod(Z ** (2 * k - 1), axis=1) * w
    s = Z.sum(axis=1) / 2; sm = np.ones_like(s); co = []
    for m in range(deg + 1): co.append(np.sum(R * sm).real / math.factorial(m)); sm = sm * s
    json.dump(dict(i0=i0, co=co, sec=time.time() - t), open(f + '.tmp', 'w')); os.replace(f + '.tmp', f)
    print('cast', i0, '%.0fs' % (time.time() - t), flush=True)
co = np.zeros(deg + 1); sec = 0
for i0 in range(N * N):
    d = json.load(open(os.path.join(CK, 'c%03d.json' % i0))); co += np.array(d['co']); sec += d['sec']
lead = cfkrs.a_k(k) * cfkrs.ks_leading(k)
out = dict(k=k, a=a, N=N, radii=radii, coeffs=co.tolist(), a_k=cfkrs.a_k(k), ks_lead=cfkrs.ks_leading(k), lead_ratio=co[-1] / lead, sec=sec, tail='v2', casti=N * N)
json.dump(out, open(OUT, 'w'), indent=1); print(k, a, N, radii, 'lead %.7f' % out['lead_ratio'], '%.0fs' % sec, flush=True)
