"""34.76: CFKRS Q_k s opraveným ocasem Eulerova součinu (cfkrs.logA tail='v2'). Argument: a (0 sudé, 1 liché).
k = 1–3: N = 24, poloměry 0,22 + 0,031 j; k = 4: N = 20, poloměry alt (0,17 … 0,245) = konečné a ref (0,22 … 0,313) = kontrola kontury."""
import sys, json, time, math; sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
a = int(sys.argv[1]); out = {}
for k, N, rad, key in [(1, 24, None, 'k1'), (2, 24, None, 'k2'), (3, 24, None, 'k3'), (4, 20, [0.17, 0.195, 0.22, 0.245], 'k4'), (4, 20, [0.22, 0.251, 0.282, 0.313], 'k4ref')]:
    t = time.time(); co = cfkrs.Qk_coeffs(k, a, N=N, radii=rad)
    out[key] = dict(coeffs=co.tolist(), N=N, radii=rad, a_k=cfkrs.a_k(k), ks_lead=cfkrs.ks_leading(k), sec=time.time() - t)
    print(k, key, co.tolist(), 'lead ratio %.7f' % (co[-1] / (cfkrs.a_k(k) * cfkrs.ks_leading(k))), '%.0fs' % (time.time() - t), flush=True)
json.dump(out, open('/home/claude/br76/data/cfkrs_Q_v2_a%d.json' % a, 'w'), indent=1)
