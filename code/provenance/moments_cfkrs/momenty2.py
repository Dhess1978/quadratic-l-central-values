"""Moments vs CFKRS with three error estimates: iid, empirical block scatter, and null models
(R10 random sample; random subsets of V with block size)."""
import json, sys, numpy as np, math
Q = json.load(open(sys.argv[1])); OUT = sys.argv[2]
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]
D = np.concatenate([v['D'] for v in V]); L = np.concatenate([v['L'] for v in V])
R = np.load('/home/claude/br76/data/R10.npz'); D100 = np.load('/home/claude/br76/data/D100.npz')
rng = np.random.default_rng(76)
def model(k, a, x):
    c = np.array(Q['k%d_a%d' % (k, a)]['coeffs']); return np.polyval(c[::-1], x)
groups = [(3, 1e4), (1e4, 3e4), (3e4, 1e5), (1e5, 3e5), (3e5, 1e6), (1e6, 3e6), (3e6, 1e7 + 1)]
res = {'groups': [], 'total': {}, 'R10': {}, 'subsets': {}, 'D100': {}}
for s, a, name in ((1, 0, 'kl'), (-1, 1, 'zap')):
    m = np.sign(D) == s; d = np.abs(D[m]); l = L[m]; x = np.log(d.astype(float))
    blk = (d - 1) // 20000
    for k in range(1, 5):
        if 'k%d_a%d' % (k, a) not in Q: continue
        mod = model(k, a, x); r = l ** k - mod
        # total
        nb = blk.max() + 1
        rb = np.bincount(blk, weights=r, minlength=nb); nbk = np.bincount(blk, minlength=nb)
        iid = math.sqrt(np.sum(np.bincount(blk, weights=r ** 2, minlength=nb) - rb ** 2 / nbk))  # sqrt(sum var*n)
        emp = math.sqrt(np.sum(rb ** 2))
        tot = dict(n=int(m.sum()), mean_data=float((l ** k).mean()), mean_cfkrs=float(mod.mean()), rel=float(r.sum() / mod.sum()),
                   err_iid_rel=iid / mod.sum(), err_emp_rel=emp / mod.sum(), z_emp=float(r.sum() / emp), z_iid=float(r.sum() / iid),
                   ratio_emp_iid=emp / iid)
        res['total']['%s_k%d' % (name, k)] = tot
        # groups (log ranges) with empirical error from block residual scatter inside group (blocks of 20000; small ranges: sub-blocks of 1000)
        for lo, hi in groups:
            g = (d >= lo) & (d < hi); sub = (d[g] - 1) // (1000 if hi <= 1e5 else 20000)
            rbg = np.bincount(sub - sub.min(), weights=r[g]); nn = int(g.sum())
            emp_g = math.sqrt(np.sum(rbg ** 2)); iid_g = float(np.std(r[g]) * math.sqrt(nn))
            res['groups'].append(dict(sign=name, k=k, lo=lo, hi=hi, n=nn, data=float((l[g] ** k).mean()), cfkrs=float(mod[g].mean()),
                                      rel=float(r[g].sum() / mod[g].sum()), rel_err_emp=emp_g / mod[g].sum(), rel_err_iid=iid_g / mod[g].sum()))
        # null model 1: R10 random sample, iid error
        mr = np.sign(R['D']) == s; xr = np.log(np.abs(R['D'][mr]).astype(float)); rr = R['L'][mr] ** k - model(k, a, xr)
        res['R10']['%s_k%d' % (name, k)] = dict(n=int(mr.sum()), rel=float(rr.sum() / model(k, a, xr).sum()), z_iid=float(rr.sum() / (np.std(rr) * math.sqrt(len(rr)))))
        # null model 2: random subsets of V (same sizes as blocks), scatter of subset residual sums vs iid
        sizes = nbk[nbk > 100]; zs = []
        for t in range(400):
            sz = int(rng.choice(sizes)); idx = rng.choice(len(r), sz, replace=False)
            zs.append(r[idx].sum() / (np.std(r) * math.sqrt(sz)))
        zs = np.array(zs)
        # contiguous blocks, same normalisation
        zb = rb[nbk > 100] / (np.std(r) * np.sqrt(nbk[nbk > 100]))
        res['subsets']['%s_k%d' % (name, k)] = dict(random_sd=float(zs.std()), blocks_sd=float(zb.std()), ratio=float(zb.std() / zs.std()))
        # D100 total
        md = np.sign(D100['D']) == s; xd = np.log(np.abs(D100['D'][md]).astype(float)); rd = D100['L'][md] ** k - model(k, a, xd)
        res['D100']['%s_k%d' % (name, k)] = dict(n=int(md.sum()), rel=float(rd.sum() / model(k, a, xd).sum()), z_iid=float(rd.sum() / (np.std(rd) * math.sqrt(len(rd)))))
        print('%s k=%d  rel %+.2e  (emp %.1e, iid %.1e, emp/iid %.3f)  z_emp %+.2f | R10 z %+.2f | blocks/random sd %.3f/%.3f | D100 rel %+.2e z %+.2f' % (
            name, k, tot['rel'], tot['err_emp_rel'], tot['err_iid_rel'], tot['ratio_emp_iid'], tot['z_emp'], res['R10']['%s_k%d' % (name, k)]['z_iid'],
            res['subsets']['%s_k%d' % (name, k)]['blocks_sd'], res['subsets']['%s_k%d' % (name, k)]['random_sd'], res['D100']['%s_k%d' % (name, k)]['rel'], res['D100']['%s_k%d' % (name, k)]['z_iid']), flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
