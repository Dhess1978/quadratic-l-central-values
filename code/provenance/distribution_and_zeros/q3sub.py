"""34.79 odv-lr-b34: polynom CFKRS Q₃ pro podrodinu D = 8d (d liché bez čtverců, kladná): lokální faktor u p = 2 jen Π_{i≤j}(1 − 2^{−1−z_i−z_j}) (χ_D(2) = 0 vždy)."""
import sys, math, json, numpy as np
sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
orig = cfkrs.logA
def logA_sub(Z, P=100000, tail='v2'):
    out = orig(Z, P, tail); k = Z.shape[1]; x = np.exp(-(0.5 + Z) * math.log(2))
    loc_full = 1 / 3 + 2 / 3 * 0.5 * (np.prod(1 / (1 - x), axis=1) + np.prod(1 / (1 + x), axis=1)); loc_sub = np.ones(Z.shape[0], complex)
    for i in range(k):
        for j in range(i, k): loc_full = loc_full * (1 - x[:, i] * x[:, j]); loc_sub = loc_sub * (1 - x[:, i] * x[:, j])
    return out - np.log(loc_full) + np.log(loc_sub)
cfkrs.logA = logA_sub
c = cfkrs.Qk_coeffs(3, 0, N=20); cfkrs.logA = orig; cf = cfkrs.Qk_coeffs(3, 0, N=20)
Q2 = json.load(open('/home/claude/br76/data/cfkrs_Q_final_v2.json'))['k3_a0']['coeffs']
print('sub', c.tolist()); print('full N=20', cf.tolist()); print('full v2', Q2)
json.dump(dict(sub=c.tolist(), full_N20=cf.tolist(), full_v2=Q2), open('/home/claude/br79/data/q3sub.json', 'w'), indent=1)
