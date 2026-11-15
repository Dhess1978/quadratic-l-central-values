"""34.85b · popisně (bez pravidla a bez verdiktu; doplněk E jmenuje ⟨L(½)⟩ proti předpovědi): momenty k = 1 … 5 hodnot L(½, χ_D) rodiny S9L (10⁹ < |D| ≤ 10⁹ + 8·10⁵)
proti úplným polynomům CFKRS — týž postup jako H0 u 10⁸ ve 34.78 (h0_h1.py): zbytek L^k − Q_k(ln|D|), chyba empiricky ze součtů po blocích šířky 2 000 v |D|,
pro srovnání jen vedoucí člen. Polynomy: pub5/repo/data/cfkrs (k = 1–4 a k = 5). Na nulách nezávisí.
Použití: python3 -I s9l_momenty.py S9L_L.npz VYSTUP.json"""
import sys, json, math, numpy as np
HV, OUT = sys.argv[1:3]; C = '/home/claude/pub5/repo/data/cfkrs/'
Q = json.load(open(C + 'cfkrs_Q_k1_k4.json'))
for a, nm in ((0, 'even'), (1, 'odd')): Q['k5_a%d' % a] = dict(coeffs=json.load(open(C + 'cfkrs_Q_k5_%s.json' % nm))['coeffs'])
def model(k, a, x, lead=False):
    c = np.array(Q['k%d_a%d' % (k, a)]['coeffs'], float)
    if lead: c = np.where(np.arange(len(c)) == len(c) - 1, c, 0.0)
    return np.polyval(c[::-1], x)
z = np.load(HV); D, L = z['D'].astype(np.int64), z['L']; out = {}; zs = []
for s, a, nm in ((1, 0, 'even'), (-1, 1, 'odd')):
    m = np.sign(D) == s; d = np.abs(D[m]); l = L[m]; x = np.log(d.astype(float)); blk = (d - 10 ** 9 - 1) // 2000
    for k in range(1, 6):
        mod = model(k, a, x); r = l ** k - mod; rb = np.bincount(blk, weights=r); emp = math.sqrt(np.sum(rb ** 2)); zk = float(r.sum() / emp)
        modl = model(k, a, x, lead=True); rl = l ** k - modl; rbl = np.bincount(blk, weights=rl)
        out['%s_k%d' % (nm, k)] = dict(n=int(m.sum()), data=float((l ** k).mean()), cfkrs=float(mod.mean()), rel=float(r.sum() / mod.sum()), rel_err=float(emp / mod.sum()), z=zk, rel_err_iid=float(r.std() * math.sqrt(len(r)) / mod.sum()),
                                       leading_term_only_rel=float(rl.sum() / modl.sum()), leading_term_only_z=float(rl.sum() / math.sqrt(np.sum(rbl ** 2))), blocks=int(blk.max() + 1))
        zs.append(abs(zk)); print('%s k=%d  data %.6g  CFKRS %.6g  rel %+.2e ± %.1e  z %+.2f | jen vedoucí člen: rel %+.3f' % (nm, k, (l ** k).mean(), mod.mean(), r.sum() / mod.sum(), emp / mod.sum(), zk, rl.sum() / modl.sum()))
out['max_abs_z'] = max(zs); out['note'] = 'descriptive, no rule recorded: same procedure as the rule for 1e8 (lab 34.78); error from sums over blocks of width 2000 in |D|'
json.dump(out, open(OUT, 'w'), indent=1); print('max |z| %.2f' % max(zs))
