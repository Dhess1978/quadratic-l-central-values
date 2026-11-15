import sys, json, time; sys.path.insert(0, '/home/claude/br76/py')
import cfkrs
out = {}
for k, N in [(1, 24), (2, 24), (3, 24), (4, int(sys.argv[1]))]:
    for a in (0, 1):
        t = time.time(); co = cfkrs.Qk_coeffs(k, a, N=N)
        out['k%d_a%d' % (k, a)] = dict(coeffs=co.tolist(), N=N, a_k=cfkrs.a_k(k), ks_lead=cfkrs.ks_leading(k), sec=time.time() - t)
        print(k, a, N, co.tolist(), '%.0fs' % (time.time() - t), flush=True)
json.dump(out, open('/home/claude/br76/data/cfkrs_Q_N%s.json' % sys.argv[1], 'w'), indent=1)
