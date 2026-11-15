"""34.82: D a nízké nuly rodin E10 (velky_2026_11_02) a E100 (e100_2026_11_05) z balíků KS_NULY → data/E_nuly.npz (D, g1, nuly v ploché řadě s ukazateli)."""
import os, gzip, json, numpy as np
SRC = ['/home/claude/br63/KS/velky_2026_11_02/E10', '/home/claude/br63/KS/e100_2026_11_05/E100']
Ds, G = [], []
for P in SRC:
    for fn in sorted(os.listdir(P)):
        if fn.endswith('.jsonl.gz'):
            for line in gzip.open(os.path.join(P, fn), 'rt'): r = json.loads(line); Ds.append(int(r['D'])); G.append([float(x) for x in r['g']])
D = np.array(Ds, np.int64); o = np.argsort(np.abs(D) * 2 + (D < 0)); D = D[o]; G = [G[i] for i in o]
assert len(np.unique(D)) == len(D), 'duplicity'
ptr = np.concatenate([[0], np.cumsum([len(g) for g in G])]); flat = np.concatenate([np.array(g) for g in G])
g1 = np.array([g[0] if len(g) else np.nan for g in G])
np.savez('/home/claude/br82/data/E_nuly.npz', D=D, g1=g1, ptr=ptr, nuly=flat)
print(len(D), 'D; sudých', int((D > 0).sum()), 'lichých', int((D < 0).sum()), '|D|', np.abs(D).min(), np.abs(D).max(), 'nul', len(flat), 'bez nul', int(np.isnan(g1).sum()))
