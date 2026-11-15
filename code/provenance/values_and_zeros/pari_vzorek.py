"""34.82: kontrola vlastní AFE proti PARI lfun na náhodném vzorku 24 D z E10 + E100 (12 sudých, 12 lichých; semínko 34822)."""
import glob, json, numpy as np, cypari2, time
pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9); pari.set_real_precision(25)
EN = np.load('/home/claude/br82/data/E_nuly.npz'); D = EN['D']; L = np.full(len(D), np.nan)
for f in glob.glob('/home/claude/br82/data/E_L2_*_2.npz'):
    z = np.load(f); L[z['i']] = z['L']
rng = np.random.default_rng(34822); i0 = rng.choice(np.nonzero(D > 0)[0], 12, replace=False); i1 = rng.choice(np.nonzero(D < 0)[0], 12, replace=False)
R = []; t0 = time.time()
for i in np.r_[i0, i1]:
    Lp = float(pari('lfun(lfuncreate(%d), 1/2)' % int(D[i]))); R.append(dict(D=int(D[i]), L_afe=float(L[i]), L_pari=Lp, rel=abs(L[i] - Lp) / Lp))
    print(R[-1], round(time.time() - t0), 's', flush=True)
json.dump(dict(vzorek=R, max_rel=max(r['rel'] for r in R)), open('/home/claude/br82/data/pari_vzorek.json', 'w'), indent=1); print('max rel', max(r['rel'] for r in R))
