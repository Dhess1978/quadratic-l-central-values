"""Moments k=1..4 by sign vs CFKRS Q_k and Keating-Snaith leading term; per-block chi^2, grouped residuals."""
import json, sys, numpy as np, math
S = json.load(open('/home/claude/br76/data/bloky_sumy.json'))
Q = json.load(open(sys.argv[1]))
groups = [[0], [1], [2, 3], list(range(4, 8)), list(range(8, 16)), list(range(16, 32)), list(range(32, 64)),
          list(range(64, 128)), list(range(128, 256)), list(range(256, 500))]
res = {}
for key, a in (('kl', 0), ('zap', 1)):
    for fam in ('V', 'D100'):
        s = S['%s_%s' % (fam, key)]; n = np.array(s['n'], float)
        X = np.array([s['x%d' % m] for m in range(11)])  # (11, nb)
        for k in range(1, 5):
            qk = Q.get('k%d_a%d' % (k, a))
            if qk is None: continue
            c = np.array(qk['coeffs']); deg = len(c) - 1
            model = (c[:, None] * X[:deg + 1]).sum(axis=0)          # sum over d of Q_k(ln|d|)
            lead = qk['a_k'] * qk['ks_lead'] * X[deg]                  # KS leading term summed
            data = np.array(s['L%d' % k]); d2 = np.array(s['L%d' % (2 * k)])
            var = d2 / n - (data / n) ** 2; sig = np.sqrt(var / n)       # error of block mean
            r = (data - model) / n
            chi2 = float(np.sum((r / sig) ** 2)); nb = len(n)
            tot = dict(mean_data=float(data.sum() / n.sum()), mean_cfkrs=float(model.sum() / n.sum()), mean_ks=float(lead.sum() / n.sum()),
                       rel_cfkrs=float((data.sum() - model.sum()) / model.sum()), rel_ks=float((data.sum() - lead.sum()) / lead.sum()),
                       z_cfkrs=float((data.sum() - model.sum()) / math.sqrt(np.sum(var * n))), chi2=chi2, nblocks=nb)
            if fam == 'V':
                gr = []
                for g in groups:
                    g = np.array(g); dn = data[g].sum(); mn = model[g].sum(); ln_ = lead[g].sum(); nn = n[g].sum()
                    err = math.sqrt(np.sum(var[g] * n[g])) / nn
                    gr.append(dict(absD_lo=1 + 20000 * int(g[0]), absD_hi=20000 * (int(g[-1]) + 1), n=int(nn), data=dn / nn, cfkrs=mn / nn, ks=ln_ / nn,
                                   rel=(dn - mn) / mn, rel_err=err / (mn / nn), z=(dn / nn - mn / nn) / err))
                tot['groups'] = gr
            res['%s_%s_k%d' % (fam, key, k)] = tot
            print('%-5s %-3s k=%d  data %.6g  CFKRS %.6g  rel %+.2e  z %+.2f  chi2/nb %.3f   KS-lead %.6g (rel %+.3f)' %
                  (fam, key, k, tot['mean_data'], tot['mean_cfkrs'], tot['rel_cfkrs'], tot['z_cfkrs'], chi2 / nb, tot['mean_ks'], tot['rel_ks']))
json.dump(res, open('/home/claude/br76/data/momenty.json', 'w'), indent=1)
