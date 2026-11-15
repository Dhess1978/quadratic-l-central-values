"""30místná kontrola mpmath: 180 náhodných + 20 nejmenších z kontrola_L_3478.json, část i z m."""
import sys, json, math, time
sys.path.insert(0, '/home/claude/br76/py')
import numpy as np, mpmath as mp, lcheck as lc
part, parts = int(sys.argv[1]), int(sys.argv[2]); mp.mp.dps = 30
rows = json.load(open('/home/claude/br78/data/kontrola_L_3478.json'))['rows']
sel = rows[0:90] + rows[1000:1090] + rows[2000:2020]
spf = lc.spf_table(200000); t0 = time.time()
import os
fn = '/home/claude/br78/data/mp3478_part%d.json' % part
out = json.load(open(fn)) if os.path.exists(fn) else []
hotovo = {r['D'] for r in out}
for i, r in enumerate(sel):
    if i % parts != part or r['D'] in hotovo: continue
    D = r['D']; q = abs(D); c = mp.mpf(1) / 4 if D > 0 else mp.mpf(3) / 4
    N = int(math.sqrt(55 * q / math.pi)) + 2; chi = lc.chi_vec(D, N, spf)
    s = mp.mpf(0); pq = mp.pi / q
    for n in np.nonzero(chi)[0]:
        n = int(n)
        if n == 0: continue
        s += int(chi[n]) * mp.gammainc(c, pq * n * n, regularized=True) / mp.sqrt(n)
    Lm = 2 * s
    out.append(dict(F=r['F'], D=D, L_mac=r['L_mac'], L_mp=mp.nstr(Lm, 25), abs_diff=float(abs(Lm - r['L_mac'])), rel_diff=float(abs(Lm - r['L_mac']) / Lm)))
    json.dump(out, open('/home/claude/br78/data/mp3478_part%d.json' % part, 'w'))
print('done', part, len(out), '%.0fs' % (time.time() - t0), flush=True)
