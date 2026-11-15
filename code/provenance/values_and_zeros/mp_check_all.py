"""30-digit mpmath check of all sampled values (kontrola_L_scipy.json rows), part i of m."""
import sys, json, math, time
sys.path.insert(0, '/home/claude/br76/py')
import numpy as np, mpmath as mp, lcheck as lc
part, parts = int(sys.argv[1]), int(sys.argv[2])
mp.mp.dps = 30
rows = json.load(open('/home/claude/br76/data/kontrola_L_scipy.json'))['rows']
spf = lc.spf_table(200000)
out = []; t0 = time.time()
for i, r in enumerate(rows):
    if i % parts != part: continue
    D = r['D']; q = abs(D); c = mp.mpf(1) / 4 if D > 0 else mp.mpf(3) / 4
    N = int(math.sqrt(55 * q / math.pi)) + 2
    chi = lc.chi_vec(D, N, spf)
    s = mp.mpf(0); pq = mp.pi / q
    for n in np.nonzero(chi)[0]:
        n = int(n)
        if n == 0: continue
        s += int(chi[n]) * mp.gammainc(c, pq * n * n, regularized=True) / mp.sqrt(n)
    Lm = 2 * s
    out.append(dict(F=r['F'], D=D, L_pari=r['L_pari'], L_mp=mp.nstr(Lm, 25), abs_diff=float(abs(Lm - r['L_pari'])), rel_diff=float(abs(Lm - r['L_pari']) / Lm)))
    if len(out) % 50 == 0:
        print(part, len(out), '%.0fs' % (time.time() - t0), flush=True)
        json.dump(out, open('/home/claude/br76/data/mp_check_part%d.json' % part, 'w'))
json.dump(out, open('/home/claude/br76/data/mp_check_part%d.json' % part, 'w'))
print('done', part, len(out), '%.0fs' % (time.time() - t0), flush=True)
