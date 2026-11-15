"""34.78 H0 odv-lq-cfkrs8 a H1 odv-lq-skal podle PRAVIDLA_3478 (pravidla zapsána před výpočtem)."""
import json, math, numpy as np
Q = json.load(open('/home/claude/br76/data/cfkrs_Q_final_v2.json'))
for a in (0, 1):
    q = json.load(open('/home/claude/br77/data/q5_a%d.json' % a)); Q['k5_a%d' % a] = dict(coeffs=q['coeffs'])
def model(k, a, x, lead=False):
    c = np.array(Q['k%d_a%d' % (k, a)]['coeffs'], float)
    if lead: c = np.where(np.arange(len(c)) == len(c) - 1, c, 0.0)
    return np.polyval(c[::-1], x)
B = np.load('/home/claude/br78/data/D100bc.npz'); DB, LB = B['D'], B['L']
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; DV = np.concatenate([v['D'] for v in V]); LV = np.concatenate([v['L'] for v in V])
out = {'H0': {}, 'H1': {}, 'H1_popis': {}}
# ---------- H0 ----------
zs = []
for s, a, nm in ((1, 0, 'kl'), (-1, 1, 'zap')):
    m = np.sign(DB) == s; d = np.abs(DB[m]); l = LB[m]; x = np.log(d.astype(float))
    blk = (d - 100100000) // 2000; blk -= blk.min()
    for k in range(1, 6):
        mod = model(k, a, x); r = l ** k - mod
        rb = np.bincount(blk, weights=r); emp = math.sqrt(np.sum(rb ** 2))
        z = float(r.sum() / emp)
        modl = model(k, a, x, lead=True); rl = l ** k - modl; rbl = np.bincount(blk, weights=rl)
        zl = float(rl.sum() / math.sqrt(np.sum(rbl ** 2)))
        out['H0']['%s_k%d' % (nm, k)] = dict(n=int(m.sum()), data=float((l ** k).mean()), cfkrs=float(mod.mean()), rel=float(r.sum() / mod.sum()),
                                          rel_err=emp / mod.sum(), z=z, nulovy_vedouci_clen_rel=float(rl.sum() / modl.sum()), nulovy_z=zl, bloku=int(blk.max() + 1))
        zs.append(abs(z))
        print('H0 %s k=%d  data %.6g  CFKRS %.6g  rel %+.2e ± %.1e  z %+.2f | jen vedoucí člen: rel %+.3f z %+.1f' % (nm, k, (l ** k).mean(), mod.mean(), r.sum() / mod.sum(), emp / mod.sum(), z, rl.sum() / modl.sum(), zl))
v0 = 'ODVOZENO' if max(zs) < 3 else ('ZAMÍTNUTO' if max(zs) >= 4 else 'NELZE ODVODIT')
out['H0']['max_abs_z'] = max(zs); out['H0']['verdikt'] = v0; print('H0 max|z| %.2f -> %s' % (max(zs), v0))
# ---------- H1 ----------
def rho(d, l, k, a, lo, hi, w):
    g = (d >= lo) & (d < hi); dd = d[g]; ll = l[g]; x = np.log(dd.astype(float))
    r = ll ** k - model(k, a, x); sig = r.std()
    b = ((dd - lo) // w).astype(np.int64); nb = int(b.max() + 1)
    S = np.bincount(b, weights=r, minlength=nb); n = np.bincount(b, minlength=nb)
    full = n > 0.5 * n.mean()  # neúplný poslední blok vynechat
    zn = S[full] / (sig * np.sqrt(n[full])); nbk = int(full.sum())
    rr = float(zn.std(ddof=1)); return rr, rr / math.sqrt(2 * (nbk - 1)), nbk
wB = [1000, 3000, 10000, 30000]; fA = math.sqrt(2e6 / 1.05e8); zz = []
for s, a, nm in ((1, 0, 'kl'), (-1, 1, 'zap')):
    mB = np.sign(DB) == s; dB = np.abs(DB[mB]); lB = LB[mB]
    mA = np.sign(DV) == s; dA = np.abs(DV[mA]); lA = LV[mA]
    for k in (1, 2):
        for w in wB:
            wa = max(1, int(round(w * fA)))
            rA, eA, nA = rho(dA, lA, k, a, 1000000, 3000001, wa)
            rB, eB, nB = rho(dB, lB, k, a, 100100000, 102000000, w)
            z = (rA - rB) / math.sqrt(eA ** 2 + eB ** 2); zz.append(abs(z))
            out['H1']['%s_k%d_w%d' % (nm, k, w)] = dict(wA=wa, rhoA=rA, errA=eA, nA=nA, rhoB=rB, errB=eB, nB=nB, z=z)
            print('H1 %s k=%d wB=%6d wA=%5d  ρA %.3f±%.3f (%d)  ρB %.3f±%.3f (%d)  z %+.2f' % (nm, k, w, wa, rA, eA, nA, rB, eB, nB, z))
        # popisně: jiné škálování a w = 1e5
        for lab, f in (('w_samo', 1.0), ('w_deleno_D', 2e6 / 1.05e8)):
            for w in wB:
                wa = max(1, int(round(w * f))); rA, eA, nA = rho(dA, lA, k, a, 1000000, 3000001, wa); rB, eB, nB = rho(dB, lB, k, a, 100100000, 102000000, w)
                out['H1_popis']['%s_%s_k%d_w%d' % (lab, nm, k, w)] = dict(wA=wa, rhoA=rA, rhoB=rB, z=(rA - rB) / math.sqrt(eA ** 2 + eB ** 2))
        rB, eB, nB = rho(dB, lB, k, a, 100100000, 102000000, 100000)
        out['H1_popis']['w1e5_%s_k%d' % (nm, k)] = dict(rhoB=rB, errB=eB, nB=nB)
v1 = 'ODVOZENO' if max(zz) < 3 else ('ZAMÍTNUTO' if max(zz) >= 4 else 'NELZE ODVODIT')
out['H1']['max_abs_z'] = max(zz); out['H1']['verdikt'] = v1; print('H1 max|z| %.2f -> %s' % (max(zz), v1))
for kk, v in out['H1_popis'].items():
    if 'z' in v: print('  popis %-28s ρA %.3f ρB %.3f z %+.2f' % (kk, v['rhoA'], v['rhoB'], v['z']))
    else: print('  popis %-28s ρB %.3f ± %.3f (%d bloků)' % (kk, v['rhoB'], v['errB'], v['nB']))
json.dump(out, open('/home/claude/br78/data/h0_h1.json', 'w'), indent=1)
