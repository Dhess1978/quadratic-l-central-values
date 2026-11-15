"""34.77: CFKRS Q_k (opravený ocas v2) pro libovolné k, a, N, poloměry. Použití: python3 qk_run.py k a N r1,r2,... VYSTUP.json"""
import sys, json, time, math; sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
k, a, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); rad = [float(x) for x in sys.argv[4].split(',')]
t = time.time(); co = cfkrs.Qk_coeffs(k, a, N=N, radii=rad)
lead = cfkrs.a_k(k) * cfkrs.ks_leading(k)
out = dict(k=k, a=a, N=N, radii=rad, coeffs=co.tolist(), a_k=cfkrs.a_k(k), ks_lead=cfkrs.ks_leading(k), lead_ratio=co[-1] / lead, sec=time.time() - t, tail='v2')
json.dump(out, open(sys.argv[5], 'w'), indent=1); print(k, a, N, rad, 'lead %.7f' % out['lead_ratio'], '%.0fs' % out['sec'], co.tolist(), flush=True)
