# Naše polynomy CFKRS (v2, k = 1–4; Q5 z konturového integrálu 34.77) proti předpovědím Aldersona a Rubinsteina, (3/π²)∫_1^X Q±(k, log t) dt
import json, mpmath as mp
mp.mp.dps = 40
Q = json.load(open('/home/claude/br76/data/cfkrs_Q_final_v2.json'))
for a in (0, 1): Q['k5_a%d' % a] = json.load(open('/home/claude/br77/data/q5_a%d.json' % a))
A = json.load(open('ar_tabulky.json'))
def I(m, L):  # ∫_0^L u^m e^u du
    return mp.e ** L * mp.fsum((-1) ** (m - j) * mp.factorial(m) / mp.factorial(j) * L ** j for j in range(m + 1)) - (-1) ** m * mp.factorial(m)
def pred(k, a, X):
    c = Q['k%d_a%d' % (k, a)]['coeffs']; L = mp.log(X)
    return 3 / mp.pi ** 2 * mp.fsum(mp.mpf(cm) * I(m, L) for m, cm in enumerate(c))
out = {}
for s, a in (('plus', 0), ('minus', 1)):
    for k in range(1, 6):
        rel = []
        for r in [r for r in A[s] if r['k'] == k]:
            p = pred(k, a, r['X']); rel.append(float(p / r['pred'] - 1))
            r['nase_pred'] = float(p); r['R_nase'] = float(r['suma'] / p)
        out['%s_k%d' % (s, k)] = rel
        print('%-5s k=%d  naše/jejich − 1: %s' % (s, k, ', '.join('%+.2e' % x for x in rel)))
json.dump(dict(rel=out, tab=A), open('srovnani_ar.json', 'w'), indent=1)
