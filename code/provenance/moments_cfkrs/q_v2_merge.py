"""34.76: sloučí cfkrs_Q_v2_a0/a1.json (opravený ocas) do cfkrs_Q_final_v2.json (k4 = poloměry alt) a q4_check_v2_a{0,1}.json (alt proti ref)."""
import json, math, numpy as np, sys; sys.path.insert(0, '/home/claude/br76/py'); import cfkrs
B = '/home/claude/br76/data/'; out = {}
for a in (0, 1):
    d = json.load(open(B + 'cfkrs_Q_v2_a%d.json' % a))
    for k in (1, 2, 3, 4):
        v = d['k%d' % k]; out['k%d_a%d' % (k, a)] = dict(coeffs=v['coeffs'], N=v['N'], a_k=v['a_k'], ks_lead=v['ks_lead'], sec=v['sec'], radii=v['radii'], tail='v2')
    ref = d['k4ref']['coeffs']; alt = d['k4']['coeffs']; lead = cfkrs.a_k(4) * cfkrs.ks_leading(4)
    res = dict(a=a, ref_radii=[0.22, 0.251, 0.282, 0.313], alt_radii=[0.17, 0.195, 0.22, 0.245], ref=ref, alt=alt, lead_ratio=[ref[-1] / lead, alt[-1] / lead], compare=[])
    for D in (1e4, 1e6, 1e7, 1e8):
        x = math.log(D); q1 = np.polyval(ref[::-1], x); q2 = np.polyval(alt[::-1], x); res['compare'].append(dict(absD=D, Q4_ref=q1, Q4_alt=q2, rel=(q2 - q1) / q1))
    json.dump(res, open(B + 'q4_check_v2_a%d.json' % a, 'w'), indent=1); print(a, res['lead_ratio'], [round(c['rel'], 7) for c in res['compare']])
json.dump(out, open(B + 'cfkrs_Q_final_v2.json', 'w'), indent=1)
old = json.load(open(B + 'cfkrs_Q_final.json'))
for key in out:
    for D in (1e4, 1e7):
        x = math.log(D); print(key, '%.0e' % D, 'v2/v1 - 1 = %.2e' % (np.polyval(out[key]['coeffs'][::-1], x) / np.polyval(old[key]['coeffs'][::-1], x) - 1))
