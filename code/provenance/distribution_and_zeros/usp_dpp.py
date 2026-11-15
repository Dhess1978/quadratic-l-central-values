"""Exact sampling of USp(2N) eigenangles theta in [0, pi] as a projection DPP with basis sqrt(2/pi) sin((n+1) theta), n < N.
Returns per sample: ln Lambda(1) = sum ln(2 - 2 cos theta_j), and x1 = theta_min * N / pi."""
import numpy as np
def sample(N, M, seed=1):
    rng = np.random.default_rng(seed)
    n = np.arange(1, N + 1)
    lnL = np.empty(M); x1 = np.empty(M)
    c = np.sqrt(2 / np.pi)
    for s in range(M):
        # orthonormal basis of the remaining span, as coefficient vectors over phi_n
        B = np.eye(N); pts = []
        for step in range(N, 0, -1):
            # density on [0,pi]: ||B^T phi(theta)||^2 / step ; bound: ||phi||^2 <= N * 2/pi
            while True:
                th = rng.uniform(0, np.pi, 64)
                ph = c * np.sin(np.outer(th, n))           # (64, N)
                v = ph @ B                                  # (64, step)
                dens = np.sum(v * v, axis=1)
                acc = rng.uniform(0, N * 2 / np.pi, 64) < dens
                if acc.any():
                    i = np.argmax(acc); t = th[i]; vt = v[i]; break
            pts.append(t)
            # remove direction vt from span B (columns): new basis orthogonal to vt within span
            if step > 1:
                w = vt / np.linalg.norm(vt)
                P = np.eye(step) - np.outer(w, w)
                U, S, _ = np.linalg.svd(P)
                B = B @ U[:, :step - 1]
        pts = np.array(pts)
        lnL[s] = np.sum(np.log(2 - 2 * np.cos(pts))); x1[s] = pts.min() * N / np.pi
    return lnL, x1
if __name__ == '__main__':
    import sys, json, time
    out = {}
    for N in range(1, 10):
        t = time.time(); lnL, x1 = sample(N, int(sys.argv[1]), seed=100 + N)
        r = np.corrcoef(lnL, x1)[0, 1]
        out[N] = dict(mean_lnL=float(lnL.mean()), var_lnL=float(lnL.var()), mean_x1=float(x1.mean()), var_x1=float(x1.var()),
                      cov=float(np.cov(lnL, x1)[0, 1]), corr=float(r), corr_log=float(np.corrcoef(lnL, np.log(x1))[0, 1]))
        print(N, out[N], '%.0fs' % (time.time() - t), flush=True)
    json.dump(out, open('/home/claude/br76/data/usp_dpp_%s.json' % sys.argv[1], 'w'), indent=1)
