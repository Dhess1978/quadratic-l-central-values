"""Link between L(1/2, chi_D) and the first zero x1 (study III families with zeros), with shuffle null model."""
import numpy as np, json, csv, math
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as f:
    rd = csv.reader(f); next(rd)
    for F, D, x1, x2, nn in rd:
        rows.setdefault(F, []).append((int(D), float(x1), float(x2)))
rng = np.random.default_rng(7601)
out = {}
for F in ['B10', 'C100', 'R10', 'D100']:
    z = np.array(rows[F]); Dz = z[:, 0].astype(np.int64); x1 = z[:, 1]; x2 = z[:, 2]
    v = np.load('/home/claude/br76/data/%s.npz' % F); order = np.argsort(v['D']); Dv = v['D'][order]; Lv = v['L'][order]
    idx = np.searchsorted(Dv, Dz); ok = (idx < len(Dv)) & (Dv[np.minimum(idx, len(Dv) - 1)] == Dz)
    print(F, 'matched', int(ok.sum()), 'of', len(Dz))
    L = Lv[idx[ok]]; Dm = Dz[ok]; sc = np.log(np.abs(Dm) / math.pi) / (2 * math.pi)
    x1 = x1[ok] * sc; x2 = x2[ok] * sc   # raw gamma -> normalized x = gamma ln(|D|/pi)/(2 pi)
    for s, name in ((1, 'kl'), (-1, 'zap')):
        m = (np.sign(Dm) == s) & np.isfinite(x1)
        y = np.log(L[m]); X = x1[m]; n = int(m.sum())
        r = float(np.corrcoef(y, X)[0, 1]); rl = float(np.corrcoef(y, np.log(X))[0, 1])
        # Spearman
        ry = np.argsort(np.argsort(y)); rx = np.argsort(np.argsort(X)); rs = float(np.corrcoef(ry, rx)[0, 1])
        shuf = [float(np.corrcoef(y, rng.permutation(X))[0, 1]) for _ in range(200)]
        # conditional mean of L in x1 quintiles
        q = np.quantile(X, [0, .2, .4, .6, .8, 1.0]); cond = []
        for i in range(5):
            mm = (X >= q[i]) & (X <= q[i + 1]); cond.append(dict(x1_lo=float(q[i]), x1_hi=float(q[i + 1]), mean_L=float(L[m][mm].mean()), mean_lnL=float(y[mm].mean())))
        out['%s_%s' % (F, name)] = dict(n=n, corr=r, corr_lnx1=rl, spearman=rs, shuffle_mean=float(np.mean(shuf)), shuffle_sd=float(np.std(shuf)),
                                        mean_x1=float(X.mean()), var_x1=float(X.var()), var_lnL=float(y.var()), cov=float(np.cov(y, X)[0, 1]),
                                        cond=cond, absD_med=float(np.median(np.abs(Dm[m]))))
        print('%-5s %s n=%6d corr(lnL,x1) %.4f  corr(lnL,ln x1) %.4f  spearman %.4f  shuffle %.4f+-%.4f  mean x1 %.4f  quintile mean L %s' % (
            F, name, n, r, rl, rs, np.mean(shuf), np.std(shuf), X.mean(), ' '.join('%.2f' % c['mean_L'] for c in cond)))
json.dump(out, open('/home/claude/br76/data/vazba_x1.json', 'w'), indent=1)
