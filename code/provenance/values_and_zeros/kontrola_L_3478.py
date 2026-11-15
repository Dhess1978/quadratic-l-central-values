# Nezávislá kontrola L(1/2, chi_D) z Macu (D100b, D100c): vlastní AFE (scipy) na 2000 náhodných + 20 nejmenších.
import sys, json, gzip, os, time, numpy as np
sys.path.insert(0, '/home/claude/br76/py'); import lcheck as lc
base = '/home/claude/br78/data_mac/STRED_D100bc_2026_11_10'
m = json.load(open(os.path.join(base, 'manifest.json')))
D_all, L_all, F_all = [], [], []
for fam, v in m['rodiny'].items():
    for b in v['bloky'].values():
        for l in gzip.open(os.path.join(base, fam, b['soubor']), 'rt'):
            r = json.loads(l); D_all.append(r['D']); L_all.append(r['L']); F_all.append(fam)
D_all = np.array(D_all); L_all = np.array(L_all); F_all = np.array(F_all)
np.savez('/home/claude/br78/data/D100bc.npz', D=D_all, L=L_all, F=F_all)
rng = np.random.default_rng(20261111)
idx = list(rng.choice(np.nonzero(F_all == 'D100b')[0], 1000, replace=False)) + list(rng.choice(np.nonzero(F_all == 'D100c')[0], 1000, replace=False))
idx += list(np.argsort(np.abs(L_all))[:20])
spf = lc.spf_table(60000); t0 = time.time(); rows = []
for i in idx:
    D = int(D_all[i]); Lp = float(L_all[i]); Lc = lc.L_half(D, spf)
    rows.append(dict(F=str(F_all[i]), D=D, L_mac=Lp, L_afe=Lc, abs=abs(Lc - Lp), rel=abs(Lc - Lp) / max(abs(Lp), 1e-300)))
rel = np.array([r['rel'] for r in rows[:2000]]); ab = np.array([r['abs'] for r in rows])
out = dict(n=len(rows), cas_s=round(time.time() - t0), max_rel_nahodne=float(rel.max()), median_rel=float(np.median(rel)), max_abs_vse=float(ab.max()),
           min_L=[(r['D'], r['L_mac'], r['L_afe']) for r in rows[2000:]], zaporne_L=int((L_all < 0).sum()), nulove_L=int((L_all == 0).sum()), rows=rows)
json.dump(out, open('/home/claude/br78/data/kontrola_L_3478.json', 'w'))
print({k: v for k, v in out.items() if k not in ('rows', 'min_L')}); print('nejmenší:', out['min_L'][:6])
