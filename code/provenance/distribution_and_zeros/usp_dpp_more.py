"""USp(2N) DPP for N = 10, 11 (extends usp_dpp_40000.json; same sampler, seeds 100 + N)."""
import sys, json, time, numpy as np; sys.path.insert(0, '/home/claude/br76/py'); from usp_dpp import sample
out = json.load(open('/home/claude/br76/data/usp_dpp_40000.json'))
for N in (10, 11):
    t = time.time(); lnL, x1 = sample(N, 40000, seed=100 + N)
    out[str(N)] = dict(mean_lnL=float(lnL.mean()), var_lnL=float(lnL.var()), mean_x1=float(x1.mean()), var_x1=float(x1.var()),
                       cov=float(np.cov(lnL, x1)[0, 1]), corr=float(np.corrcoef(lnL, x1)[0, 1]), corr_log=float(np.corrcoef(lnL, np.log(x1))[0, 1]))
    print(N, out[str(N)], '%.0fs' % (time.time() - t), flush=True)
json.dump(out, open('/home/claude/br76/data/usp_dpp_40000_N11.json', 'w'), indent=1)
