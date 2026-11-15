"""34.79 odv-lr-hough podle PRAVIDLA_3479: P(L < 0,03 | x₁ < 0,2) / P(L < 0,03) v 8 skupinách proti USp(2N_x1).
USp: Haarovy matice USp(2N) (kvaternionový Gram–Schmidt), N = 2, 3, 4 a 5 po 400 000 vzorcích, lineárně mezi celými N; práh pro Λ(1) volen tak, aby
P(Λ < ε′) = P_data(L < 0,03) (shoda kvantilu odstraní aritmetické měřítko). Odchylka od pravidla: místo „přesně (Jacobiho soubor)“ simulace s chybou v z."""
import json, math, csv, numpy as np, scipy.special as sp
rng = np.random.default_rng(34793)
def haar_usp(N, S):
    n2 = 2 * N; J = np.zeros((n2, n2)); J[:N, N:] = np.eye(N); J[N:, :N] = -np.eye(N); C = []
    for j in range(N):
        v = (rng.standard_normal((S, n2)) + 1j * rng.standard_normal((S, n2))) / math.sqrt(2)
        for u in C: v = v - np.sum(np.conj(u) * v, 1)[:, None] * u
        v = v / np.linalg.norm(v, axis=1)[:, None]; C += [v, -(np.conj(v) @ J.T)]
    U = np.stack(C[0::2] + C[1::2], 2); th = np.angle(np.linalg.eigvals(U)); th = np.sort(np.where(th > 0, th, 9.0), 1)[:, :N]
    return np.sum(np.log(2 - 2 * np.cos(th)), 1), th[:, 0] * N / math.pi
SIM = {}
for N in (2, 3, 4, 5):
    lnL, x1 = [], []
    for _ in range(4):
        a, b = haar_usp(N, 100000); lnL.append(a); x1.append(b)
    SIM[N] = (np.concatenate(lnL), np.concatenate(x1)); print('USp N', N, 'P(x1<0,2)', round(float(np.mean(SIM[N][1] < 0.2)), 5), flush=True)
def usp_ratio(N, q):
    lnL, x1 = SIM[N]; thr = np.quantile(lnL, q); a = lnL < thr; c = x1 < 0.2; pa = a.mean(); pc = c.mean(); pac = (a & c).mean()
    r = (pac / pc) / pa; n = len(lnL); e = r * math.sqrt(max(1 / max((a & c).sum(), 1) + 1 / c.sum(), 0)); return r, e
H2 = json.load(open('/home/claude/br78/data/h2_vaha.json'))
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as fh:
    rd = csv.reader(fh); next(rd)
    for F, d0, g1, g2, nn in rd: rows.setdefault(F, []).append((int(d0), float(g1)))
out = {'skupiny': {}}; zs = []
for F in ('B10', 'C100', 'R10', 'D100'):
    f = np.load('/home/claude/br76/data/' + F + '.npz'); lv = dict(zip(f['D'].tolist(), f['L'].tolist()))
    for nm, sg, a in (('kl', 1, 0), ('zap', -1, 1)):
        pr = [(d0, g1, lv[d0]) for d0, g1 in rows[F] if (d0 > 0) == (sg > 0) and d0 in lv and math.isfinite(g1)]
        Dv = np.array([p[0] for p in pr], np.int64); t = np.array([p[1] for p in pr]); L = np.array([p[2] for p in pr])
        x = ((t / 2) * np.log(np.abs(Dv) / math.pi) + np.imag(sp.loggamma((0.5 + a + 1j * t) / 2))) / math.pi
        sm = L < 0.03; lo = x < 0.2; q = float(sm.mean()); pc = float(lo.mean()); nac = int((sm & lo).sum())
        R = (nac / lo.sum()) / q; eR = R * math.sqrt(1 / max(nac, 1) + 1 / lo.sum() + 1 / sm.sum())
        Nx = H2['skupiny']['%s_%s' % (F, nm)]['N_x1']; i = min(max(int(math.floor(Nx)), 2), 4); tt = Nx - i
        r0, e0 = usp_ratio(i, q); r1, e1 = usp_ratio(i + 1, q); Ru = r0 * (1 - tt) + r1 * tt; eu = math.hypot(e0 * (1 - tt), e1 * tt)
        Lp = rng.permutation(L); Rn = (np.sum((Lp < 0.03) & lo) / lo.sum()) / q
        z = (R - Ru) / math.hypot(eR, eu); zs.append(z)
        out['skupiny']['%s_%s' % (F, nm)] = dict(n=len(L), q=q, p_nizka=pc, n_oboje=nac, R=R, eR=eR, N_x1=Nx, R_usp=Ru, e_usp=eu, z=z, R_nul=float(Rn))
        print('%-9s n %6d P(L<0,03) %.4f P(x₁<0,2) %.4f oboje %4d  R %.2f ± %.2f  USp(2·%.2f) %.2f ± %.2f  z %+.1f  promíchání %.2f' % ('%s_%s' % (F, nm), len(L), q, pc, nac, R, eR, Nx, Ru, eu, z, Rn), flush=True)
zm = max(abs(z) for z in zs); out['zmax'] = zm; out['verdikt'] = 'O' if zm < 3 else ('Z' if zm >= 4 else 'N'); print('max|z| %.1f verdikt %s' % (zm, out['verdikt']))
json.dump(out, open('/home/claude/br79/data/hough.json', 'w'), indent=1)
