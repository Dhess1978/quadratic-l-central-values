"""Q_4 robustness: recompute with different radii (same N=20) and compare Q_4(ln|D|) at |D| = 1e4, 1e6, 1e7, 1e8."""
import sys, json, time, math, numpy as np; sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
a = int(sys.argv[1])
ref = json.load(open('/home/claude/br76/data/cfkrs_Q_N20.json'))['k4_a%d' % a]['coeffs']
t = time.time()
alt = cfkrs.Qk_coeffs(4, a, N=20, radii=[0.17, 0.195, 0.22, 0.245]).tolist()
res = {'a': a, 'ref_radii': [0.22, 0.251, 0.282, 0.313], 'alt_radii': [0.17, 0.195, 0.22, 0.245], 'ref': ref, 'alt': alt, 'sec': time.time() - t, 'compare': []}
for D in (1e4, 1e6, 1e7, 1e8):
    x = math.log(D); q1 = np.polyval(ref[::-1], x); q2 = np.polyval(alt[::-1], x)
    res['compare'].append(dict(absD=D, Q4_ref=q1, Q4_alt=q2, rel=(q2 - q1) / q1))
    print(a, '%.0e' % D, q1, q2, '%.2e' % ((q2 - q1) / q1), flush=True)
lead = cfkrs.a_k(4) * cfkrs.ks_leading(4)
print('lead ratio ref %.6f alt %.6f' % (ref[-1] / lead, alt[-1] / lead))
res['lead_ratio'] = [ref[-1] / lead, alt[-1] / lead]
json.dump(res, open('/home/claude/br76/data/q4_check_a%d.json' % a, 'w'), indent=1)
