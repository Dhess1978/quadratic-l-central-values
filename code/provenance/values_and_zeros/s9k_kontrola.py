"""34.85b · kontrola L(½, χ_D) u 10⁹ podle PRAVIDLA_3485.md (část C a doplněk E, zapsáno před výpočtem):
PARI lfun z Macu (rodina S9K = 2 000 náhodných D z S9L, ks_nuly_3485.py --hodnota-stred) proti vlastní funkcionální rovnici (afe_L3.py, data/S9L_L_*_2.npz).
Hodnoty z cloudu se použijí, když největší relativní rozdíl na všech D vzorku je pod 1e-9 (a každé D vzorku leží v S9L).
Použití: python3 -I s9k_kontrola.py SLOZKA_S9K S9L_D.npz PREFIX_HODNOT VYSTUP.json [S9L_L.npz]   (SLOZKA = …/S9K_hodnoty_…/S9K)"""
import sys, os, gzip, json, numpy as np
P, DN, PREF, OUT = sys.argv[1:5]
pari = {}
for fn in sorted(os.listdir(P)):
    if fn.endswith('.jsonl.gz'):
        for line in gzip.open(os.path.join(P, fn), 'rt'):
            r = json.loads(line); d = int(r['D']); assert d not in pari; pari[d] = float(r['L'])
DE = np.load(DN)['D'].astype(np.int64); afe = {}; casti = sorted(f for f in os.listdir(os.path.dirname(PREF) or '.') if f.startswith(os.path.basename(PREF) + '_') and f.endswith('.npz'))
for f in casti:
    z = np.load(os.path.join(os.path.dirname(PREF) or '.', f))
    for i, L in zip(z['i'].tolist(), z['L'].tolist()): afe[int(DE[i])] = L
spol = sorted(set(pari) & set(afe)); jen_pari = sorted(set(pari) - set(afe))
Lp = np.array([pari[d] for d in spol]); La = np.array([afe[d] for d in spol]); rel = np.abs(Lp - La) / np.maximum(np.abs(Lp), 1e-300); i = int(np.argmax(rel)); kl = np.array(spol) > 0
out = dict(n_pari=len(pari), n_afe=len(afe), n_S9L=int(len(DE)), n_spolecne=len(spol), n_jen_pari=len(jen_pari), absmin=int(min(abs(d) for d in pari)), absmax=int(max(abs(d) for d in pari)),
           max_rel=float(rel.max()), D_max_rel=int(spol[i]), L_pari_max=float(Lp[i]), L_afe_max=float(La[i]), median_rel=float(np.median(rel)), max_abs=float(np.abs(Lp - La).max()),
           n_kl=int(kl.sum()), n_zap=int((~kl).sum()), L_min_pari=float(Lp.min()), L_min_afe=float(min(afe.values())), afe_nekladne=int(sum(v <= 0 for v in afe.values())),
           kvantily_rel=[float(q) for q in np.quantile(rel, [0.5, 0.9, 0.99, 0.999])], afe_uplne=bool(len(afe) == len(DE)))
out['kontrola'] = 'O' if (out['max_rel'] < 1e-9 and out['n_jen_pari'] == 0 and out['afe_uplne']) else 'X'
print(json.dumps(out, ensure_ascii=False)); json.dump(out, open(OUT, 'w'), indent=1)
if out['kontrola'] == 'O' and len(sys.argv) > 5:
    L = np.array([afe[int(d)] for d in DE]); np.savez(sys.argv[5], D=DE, L=L); print('hodnoty S9L uloženy', sys.argv[5], len(L))
