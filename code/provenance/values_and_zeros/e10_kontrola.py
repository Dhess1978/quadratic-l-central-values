"""34.84 · kontrola L(½, χ_D) u 10¹⁰ podle PRAVIDLA_3484.md (oddíl 2, zapsáno před výpočtem):
PARI lfun z Macu (STRED_E10_2026_11_12, ks_nuly_3482.py --hodnota-stred) proti vlastní funkcionální rovnici 34.82 (afe_L2.py, br82/data/E_L2_*_2.npz).
Kontrola v pořádku, když největší relativní rozdíl < 1e-9 u všech společných D.
Použití: python3 -I e10_kontrola.py SLOZKA_E10 VYSTUP.json   (SLOZKA = …/STRED_E10_2026_11_12/E10)"""
import sys, os, gzip, json, numpy as np
P, OUT = sys.argv[1], sys.argv[2]
pari = {}
for fn in sorted(os.listdir(P)):
    if fn.endswith('.jsonl.gz'):
        for line in gzip.open(os.path.join(P, fn), 'rt'):
            r = json.loads(line); d = int(r['D']); assert d not in pari; pari[d] = float(r['L'])
E = np.load('/home/claude/br82/data/E_nuly.npz'); DE = E['D'].astype(np.int64)
afe = {}
for a in (0, 1):
    z = np.load('/home/claude/br82/data/E_L2_%d_2.npz' % a)
    for i, L in zip(z['i'].tolist(), z['L'].tolist()): afe[int(DE[i])] = L
spol = sorted(set(pari) & set(afe)); jen_pari = sorted(set(pari) - set(afe))
Lp = np.array([pari[d] for d in spol]); La = np.array([afe[d] for d in spol])
rel = np.abs(Lp - La) / np.maximum(np.abs(Lp), 1e-300); ab = np.abs(Lp - La)
i = int(np.argmax(rel)); kl = np.array(spol) > 0
out = dict(n_pari=len(pari), n_spolecne=len(spol), n_jen_pari=len(jen_pari), absmin=int(min(abs(d) for d in pari)), absmax=int(max(abs(d) for d in pari)),
           max_rel=float(rel.max()), D_max_rel=int(spol[i]), L_pari_max=float(Lp[i]), L_afe_max=float(La[i]), median_rel=float(np.median(rel)), max_abs=float(ab.max()),
           max_rel_kl=float(rel[kl].max()), max_rel_zap=float(rel[~kl].max()), n_kl=int(kl.sum()), n_zap=int((~kl).sum()),
           L_min=float(Lp.min()), pod_1e6=int((Lp < 1e-6).sum()), zaporne=int((Lp < 0).sum()), kvantily_rel=[float(q) for q in np.quantile(rel, [0.5, 0.9, 0.99, 0.999])])
out['kontrola'] = 'O' if (out['max_rel'] < 1e-9 and out['n_jen_pari'] == 0) else 'X'
print(json.dumps({k: v for k, v in out.items()}, ensure_ascii=False))
json.dump(out, open(OUT, 'w'), indent=1)
