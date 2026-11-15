"""34.77 STŘED V · DOPLŇKY — souhrnné součty a měření pro list (vstupy: br76/data/*.npz, nuly_sklad, CFKRS v2 + Q₅). Výstup data/an3477.json.
Pravidla: br77/PRAVIDLA_3477.md (zapsána před výpočtem)."""
import json, math, csv, sys, numpy as np
from scipy.special import loggamma
B6 = '/home/claude/br76/data/'; B7 = '/home/claude/br77/data/'
V = [np.load(B6 + '%s.npz' % f) for f in ('V5', 'V10')]; D = np.concatenate([v['D'] for v in V]).astype(np.int64); L = np.concatenate([v['L'] for v in V])
FAM = {f: np.load(B6 + '%s.npz' % f) for f in ('B10', 'C100', 'R10', 'D100')}
Q = json.load(open(B6 + 'cfkrs_Q_final_v2.json'))
import os
MAM_Q5 = all(os.path.exists(B7 + 'q5_a%d.json' % a) for a in (0, 1))
for a in ((0, 1) if MAM_Q5 else ()):
    q5 = json.load(open(B7 + 'q5_a%d.json' % a)); Q['k5_a%d' % a] = dict(coeffs=q5['coeffs'], N=q5['N'], a_k=q5['a_k'], ks_lead=q5['ks_lead'], lead_ratio=q5['lead_ratio'], radii=q5['radii'])
SG = (('kl', 1, 0), ('zap', -1, 1)); out = {}
def model(k, a, x): c = np.array(Q['k%d_a%d' % (k, a)]['coeffs']); return np.polyval(c[::-1], x)
# ---- 1) bloky (n, L^1..L^10, x^1..x^15, min) a součty pro přesnou iid chybu (X^0..30, L^k x^j pro k ≤ 5, j ≤ 15)
out['bl'] = {}; out['sum'] = {}
def sums(dd, ll):
    x = np.log(np.abs(dd).astype(float)); return dict(n=int(len(dd)), L=[float(np.sum(ll ** k)) for k in range(1, 11)], X=[float(np.sum(x ** m)) for m in range(31)], LX=[[float(np.sum(ll ** k * x ** j)) for j in range(16)] for k in range(1, 6)])
for nm, s, a in SG:
    m = np.sign(D) == s; d = np.abs(D[m]); l = L[m]; x = np.log(d.astype(float)); b = (d - 1) // 20000; nb = int(b.max()) + 1
    cols = [np.bincount(b, minlength=nb).astype(float)] + [np.bincount(b, weights=l ** k, minlength=nb) for k in range(1, 11)] + [np.bincount(b, weights=x ** j, minlength=nb) for j in range(1, 16)]
    mn = np.full(nb, np.inf); np.minimum.at(mn, b, l); cols.append(mn)
    out['bl'][nm] = np.column_stack(cols).tolist(); out['sum']['V_' + nm] = sums(D[m], L[m])
    for F in ('R10', 'D100'):
        f = FAM[F]; mm = np.sign(f['D']) == s; out['sum'][F + '_' + nm] = sums(f['D'][mm], f['L'][mm])
# ---- 2) k = 5: celkem, pásma, nulové modely (jako momenty2.py)
groups = [(3, 1e4), (1e4, 3e4), (3e4, 1e5), (1e5, 3e5), (3e5, 1e6), (1e6, 3e6), (3e6, 1e7 + 1)]
out['k5'] = {'total': {}, 'groups': [], 'R10': {}, 'D100': {}, 'subsets': {}}; rng = np.random.default_rng(77)
for nm, s, a in (SG if MAM_Q5 else ()):
    m = np.sign(D) == s; d = np.abs(D[m]); l = L[m]; x = np.log(d.astype(float)); blk = (d - 1) // 20000; mod = model(5, a, x); r = l ** 5 - mod
    nb = blk.max() + 1; rb = np.bincount(blk, weights=r, minlength=nb); nbk = np.bincount(blk, minlength=nb)
    iid = math.sqrt(np.sum(np.bincount(blk, weights=r ** 2, minlength=nb) - rb ** 2 / nbk)); emp = math.sqrt(np.sum(rb ** 2))
    out['k5']['total'][nm] = dict(n=int(m.sum()), mean_data=float((l ** 5).mean()), mean_cfkrs=float(mod.mean()), rel=float(r.sum() / mod.sum()), err_emp_rel=emp / mod.sum(), err_iid_rel=iid / mod.sum(), z_emp=float(r.sum() / emp), z_iid=float(r.sum() / iid))
    for lo, hi in groups:
        g = (d >= lo) & (d < hi); sub = (d[g] - 1) // (1000 if hi <= 1e5 else 20000); rbg = np.bincount(sub - sub.min(), weights=r[g])
        out['k5']['groups'].append(dict(s=nm, lo=lo, hi=hi, n=int(g.sum()), data=float((l[g] ** 5).mean()), cfkrs=float(mod[g].mean()), rel=float(r[g].sum() / mod[g].sum()), e=math.sqrt(np.sum(rbg ** 2)) / mod[g].sum(), ei=float(np.std(r[g]) * math.sqrt(g.sum())) / mod[g].sum()))
    for F in ('R10', 'D100'):
        f = FAM[F]; mm = np.sign(f['D']) == s; xr = np.log(np.abs(f['D'][mm]).astype(float)); rr = f['L'][mm] ** 5 - model(5, a, xr)
        out['k5'][F][nm] = dict(n=int(mm.sum()), rel=float(rr.sum() / model(5, a, xr).sum()), z_iid=float(rr.sum() / (np.std(rr) * math.sqrt(len(rr)))))
    sizes = nbk[nbk > 100]; zs = [r[rng.choice(len(r), int(sz), replace=False)].sum() / (np.std(r) * math.sqrt(sz)) for sz in rng.choice(sizes, 400)]
    zb = rb[nbk > 100] / (np.std(r) * np.sqrt(nbk[nbk > 100])); out['k5']['subsets'][nm] = dict(random_sd=float(np.std(zs)), blocks_sd=float(zb.std()), ratio=float(zb.std() / np.std(zs)))
    lead = Q['k5_a%d' % a]['a_k'] * Q['k5_a%d' % a]['ks_lead']; out['k5']['total'][nm]['ks'] = float((l ** 5).sum() / (lead * np.sum(x ** 15)) - 1)
print('k5', {k: (round(v['rel'], 6), round(v['z_emp'], 2)) for k, v in out['k5']['total'].items()}, flush=True)
# ---- 3) první nula: rozvinutí θ_D(γ)/π, třídy χ_D(2), χ_D(3), χ_D(5)
def chi_p(Dv, p):
    if p == 2: r = Dv % 8; return np.where(Dv % 2 == 0, 0, np.where((r == 1) | (r == 7), 1, -1))
    r = Dv % p; e = np.array([pow(int(t), (p - 1) // 2, p) for t in range(p)]); v = e[r]; return np.where(r == 0, 0, np.where(v == 1, 1, -1))
rows = {}
with open('/home/claude/br76/nuly_sklad/prvni_nuly_B10_C100_D100_R10.csv') as fh:
    rd = csv.reader(fh); next(rd)
    for F, d0, g1, g2, nn in rd: rows.setdefault(F, []).append((int(d0), float(g1)))
out['x1'] = {}
for F in ('B10', 'C100', 'R10', 'D100'):
    f = FAM[F]; o = np.argsort(f['D']); Dv = f['D'][o]; Lv = f['L'][o]; z = np.array(rows[F]); Dz = z[:, 0].astype(np.int64); g1 = z[:, 1]
    idx = np.searchsorted(Dv, Dz); ok = (idx < len(Dv)) & (Dv[np.minimum(idx, len(Dv) - 1)] == Dz) & np.isfinite(g1); Lm = Lv[idx[ok]]; Dm = Dz[ok]; t = g1[ok]
    for nm, s, a in SG:
        m = np.sign(Dm) == s; d = np.abs(Dm[m]).astype(float); tt = t[m]; y = np.log(Lm[m])
        x = ((tt / 2) * np.log(d / math.pi) + np.imag(loggamma((0.5 + a + 1j * tt) / 2))) / math.pi
        c2, c3, c5 = chi_p(Dm[m], 2), chi_p(Dm[m], 3), chi_p(Dm[m], 5); st = (c2 + 1) * 9 + (c3 + 1) * 3 + (c5 + 1)
        S = []
        for k in range(27):
            mk = st == k
            if mk.sum() < 2: continue
            S.append([k, int(mk.sum()), float(x[mk].sum()), float(y[mk].sum()), float((x[mk] ** 2).sum()), float((y[mk] ** 2).sum()), float((x[mk] * y[mk]).sum())])
        xr = x.copy(); yr = y.copy()
        for k in range(27):
            mk = st == k
            if mk.any(): xr[mk] -= x[mk].mean(); yr[mk] -= y[mk].mean()
        out['x1']['%s_%s' % (F, nm)] = dict(n=int(m.sum()), Sx=float(x.sum()), Sy=float(y.sum()), Sxx=float((x * x).sum()), Syy=float((y * y).sum()), Sxy=float((x * y).sum()), Slog=float(np.sum(np.log(d))),
                                           corr=float(np.corrcoef(x, y)[0, 1]), corr_tr=float(np.sum(xr * yr) / math.sqrt(np.sum(xr * xr) * np.sum(yr * yr))), tridy=S, mx1=float(x.mean()), vx1=float(x.var()))
    print(F, {k: (round(v['corr'], 4), round(v['corr_tr'], 4), round(v['mx1'], 4)) for k, v in out['x1'].items() if k.startswith(F)}, flush=True)
# ---- 4) ln L po pásmech: mocninné součty, zlomkové momenty, chyby κ₃ (iid a jackknife po blocích)
out['lnL'] = []
for lo, hi in [(1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 3e6), (3e6, 1e7 + 1)]:
    for nm, s, a in SG:
        m = (np.sign(D) == s) & (np.abs(D) >= lo) & (np.abs(D) < hi); y = np.log(L[m]); d = np.abs(D[m]); n = len(y)
        def k3(v): c = v - v.mean(); return float(np.mean(c ** 3))
        c = y - y.mean(); m2 = np.mean(c ** 2); m3 = np.mean(c ** 3); m4 = np.mean(c ** 4); m6 = np.mean(c ** 6)
        e3_iid = math.sqrt(max(m6 - m3 ** 2 - 6 * m4 * m2 + 9 * m2 ** 3, 0) / n)   # delta metoda pro centrální 3. moment
        bb = (d - 1) // 20000; ub = np.unique(bb); e3_jk = None
        if len(ub) >= 10:
            jk = np.array([k3(y[bb != t]) for t in ub]); e3_jk = math.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2))
        out['lnL'].append(dict(s=nm, lo=lo, hi=hi, n=n, Sy=[float(np.sum(y ** p)) for p in range(1, 5)], Sx=float(np.sum(np.log(d.astype(float)))),
                               Lh=float(np.sum(L[m] ** 0.5)), Lmh=float(np.sum(L[m] ** -0.5)), Lh2=float(np.sum(L[m])), Lmh2=float(np.sum(L[m] ** -1.0)), e3_iid=e3_iid, e3_jk=e3_jk))
print('lnL hotovo', flush=True)
# ---- 5) malé hodnoty po třídách χ_D(2), χ_D(3) (|D| ≥ 10⁶)
out['male_tr'] = []
for lo, hi in [(1e6, 3e6), (3e6, 1e7 + 1)]:
    for nm, s, a in SG:
        m = (np.sign(D) == s) & (np.abs(D) >= lo) & (np.abs(D) < hi); Dm = D[m]; c2, c3 = chi_p(Dm, 2), chi_p(Dm, 3)
        for v2 in (-1, 0, 1):
            for v3 in (-1, 0, 1):
                mk = (c2 == v2) & (c3 == v3); out['male_tr'].append(dict(s=nm, lo=lo, hi=hi, c2=v2, c3=v3, n=int(mk.sum()), obs=int((L[m][mk] < 0.03).sum()), Sx=float(np.sum(np.log(np.abs(Dm[mk]).astype(float))))))
print('male_tr', [(t['s'], t['c2'], t['c3'], t['n'], t['obs']) for t in out['male_tr'] if t['lo'] == 3e6], flush=True)
out['Q5'] = {a: Q['k5_a%d' % a] for a in (0, 1)} if MAM_Q5 else None
json.dump(out, open(B7 + 'an3477.json', 'w')); print('ok')
