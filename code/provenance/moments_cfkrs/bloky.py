"""Per-block power sums for study V: V5+V10 (all fundamental D, 3<=|D|<=1e7) in blocks of 20000 in |D|; D100 separately."""
import numpy as np, json
D = np.concatenate([np.load('/home/claude/br76/data/V5.npz')['D'], np.load('/home/claude/br76/data/V10.npz')['D']])
L = np.concatenate([np.load('/home/claude/br76/data/V5.npz')['L'], np.load('/home/claude/br76/data/V10.npz')['L']])
out = {}
def sums(Dsel, Lsel, blk):
    x = np.log(np.abs(Dsel).astype(float)); nb = blk.max() + 1
    r = {'n': np.bincount(blk, minlength=nb).tolist()}
    for k in range(1, 9): r['L%d' % k] = np.bincount(blk, weights=Lsel ** k, minlength=nb).tolist()
    for m in range(0, 11): r['x%d' % m] = np.bincount(blk, weights=x ** m, minlength=nb).tolist()
    r['absD_mid'] = [float(np.median(np.abs(Dsel[blk == b]))) if (blk == b).any() else 0.0 for b in range(nb)]
    return r
blk_all = (np.abs(D) - 1) // 20000
for s, name in ((1, 'kl'), (-1, 'zap')):
    m = np.sign(D) == s
    out['V_' + name] = sums(D[m], L[m], blk_all[m])
D1 = np.load('/home/claude/br76/data/D100.npz'); Dd, Ld = D1['D'], D1['L']
blk = (np.abs(Dd) - Dd.__abs__().min()) // 2000   # 50 blocks of 2000 in |D|
for s, name in ((1, 'kl'), (-1, 'zap')):
    m = np.sign(Dd) == s
    out['D100_' + name] = sums(Dd[m], Ld[m], blk[m])
json.dump(out, open('/home/claude/br76/data/bloky_sumy.json', 'w'))
for k, v in out.items(): print(k, len(v['n']), sum(v['n']), 'mean L %.5f' % (sum(v['L1']) / sum(v['n'])))
# duplicates consistency: B10, C100, R10 vs V5/V10
ref = dict(zip(D.tolist(), L.tolist()))
for F in ['B10', 'C100', 'R10']:
    f = np.load('/home/claude/br76/data/%s.npz' % F); dd = [abs(ref[d] - l) / l for d, l in zip(f['D'].tolist(), f['L'].tolist()) if d in ref]
    print(F, 'overlap', len(dd), 'of', len(f['D']), 'max rel diff %.2e' % max(dd))
