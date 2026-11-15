#!/usr/bin/env python3
"""Independent check of study V. Shares no code with code/provenance: reads the published central values (data/central_values)
and the results file, and recomputes from scratch
  1. the set of fundamental discriminants 3 <= |D| <= 10^7 (own sieve) against the published values,
  2. the moments k = 1..5 against the CFKRS polynomials by band and in total, and at 10^8,
  3. the parity difference of the first moments against A0*pi/2,
  4. the cancellation of consecutive blocks of 20,000 against random subsets (own seed),
  5. the maximum-likelihood exponent of small values (eps0 = 0.03) and the smallest values,
  6. the checks of the values at 10^10 and at 10^9 (own functional equation against PARI; at 10^9 also the set of discriminants by an own segmented sieve and, descriptively, the moments against CFKRS),
  7. the chi2 and amplitudes of the corrected weighted density, the sizes of the families S7 and S9L (own sieves) and the sums of the explicit formula on C100 and S7 (Section 4.6).
usage: python3 -I independent_check_v.py REPOSITORY_ROOT"""
import sys, os, json, math, numpy as np
from scipy import optimize, stats, special
ROOT = sys.argv[1]; R = json.load(open(os.path.join(ROOT, 'data/results/study_v_results.json'))); CV = os.path.join(ROOT, 'data/central_values')
def load(name): z = np.load(os.path.join(CV, name + '.npz')); return z['D'], z['L']
ok = []
def check(name, cond, info): ok.append(bool(cond)); print(('PASS ' if cond else 'FAIL ') + name + ': ' + info, flush=True)
# ---------------------------------------------------------------- 1. sieve of fundamental discriminants
parts = [p for p in sorted(os.listdir(CV)) if p.startswith('V_part')]
D = np.concatenate([load(p[:-4])[0] for p in parts]); L = np.concatenate([load(p[:-4])[1] for p in parts])
X = 10 ** 7; sf = np.ones(X + 1, bool); sf[0] = False
for p in range(2, int(X ** 0.5) + 1): sf[p * p::p * p] = False
n = np.arange(X + 1); fund = set()
# D = 1 mod 4 squarefree, or D = 4m with m = 2, 3 mod 4 squarefree (both signs)
for sgn in (1, -1):
    m = n[sf & (n >= 1)]; d = sgn * m; ok1 = d[(d % 4 == 1) & (np.abs(d) >= 3)]
    m2 = n[sf & (n >= 1) & (4 * n <= X)]; d2 = sgn * m2; ok2 = 4 * d2[(d2 % 4 == 2) | (d2 % 4 == 3)]
    fund |= set(ok1.tolist()) | set(ok2.tolist())
check('sieve', len(fund) == len(D) == len(set(D.tolist())) and fund == set(D.tolist()), 'own sieve %d fundamental discriminants, published values %d' % (len(fund), len(D)))
check('positive values', L.min() > 0, 'smallest %.6g at D = %d' % (L.min(), D[np.argmin(L)]))
# ---------------------------------------------------------------- 2. moments against CFKRS
P = R['cfkrs_polynomials']
def model(k, s, x): c = np.array(P['k%d_%s' % (k, s)]); return np.polyval(c[::-1], x)
BANDS = [(3, 1e4), (1e4, 3e4), (3e4, 1e5), (1e5, 3e5), (3e5, 1e6), (1e6, 3e6), (3e6, 1e7 + 1)]
pub = {(g['s'], g['k'], g['lo']): g for g in R['moments']['bands']}
dmax = 0.0; tot = {}
for s, sg in (('even', 1), ('odd', -1)):
    m = np.sign(D) == sg; a = np.abs(D[m]).astype(float); l = L[m]; x = np.log(a)
    for k in range(1, 6):
        mod = model(k, s, x); tot[(s, k)] = (l ** k).sum() / mod.sum() - 1
        dmax = max(dmax, abs(tot[(s, k)] - R['moments']['total']['%s_k%d' % (s, k)]['rel']))
        for lo, hi in BANDS:
            g = (a >= lo) & (a < hi); rel = (l[g] ** k).sum() / mod[g].sum() - 1; dmax = max(dmax, abs(rel - pub[(s, k, lo)]['rel']))
check('moments', dmax < 1e-9, 'k = 1..5, both signs, 7 bands and totals: relative differences = published to %.1e; totals %s' % (dmax, ', '.join('%s k%d %+.1e' % (s, k, v) for (s, k), v in sorted(tot.items()))))
Db, Lb = np.concatenate([load('D100b')[0], load('D100c')[0]]), np.concatenate([load('D100b')[1], load('D100c')[1]]); d8 = 0.0
for s, sg in (('even', 1), ('odd', -1)):
    m = np.sign(Db) == sg; x = np.log(np.abs(Db[m]).astype(float))
    for k in range(1, 6): d8 = max(d8, abs((Lb[m] ** k).sum() / model(k, s, x).sum() - 1 - R['moments']['at_1e8']['%s_k%d' % (s, k)]['rel']))
check('moments at 1e8', d8 < 1e-9, 'D100b + D100c (%d values): relative differences = published to %.1e' % (len(Db), d8))
# ---------------------------------------------------------------- 3. parity
dif = L[D < 0].mean() - L[D > 0].mean(); A0 = R['moments']['k1_constants']['A0']; pred = A0 * math.pi / 2
check('parity', abs(dif - R['moments']['parity']['dif']) < 1e-12 and abs(pred - R['moments']['parity']['pred']) < 1e-6, 'mean odd - mean even = %.5f, A0*pi/2 = %.5f' % (dif, pred))
# ---------------------------------------------------------------- 4. consecutive blocks against random subsets (k = 1)
rng = np.random.default_rng(20261113); rat = {}
for s, sg in (('even', 1), ('odd', -1)):
    m = np.sign(D) == sg; a = np.abs(D[m]); l = L[m]; r = l - model(1, s, np.log(a.astype(float)))
    blk = (a - 1) // 20000; sums = np.bincount(blk, weights=r); cnt = np.bincount(blk); use = cnt > 100
    sd_blocks = np.sqrt(np.mean(sums[use] ** 2 / cnt[use]))
    rs = []
    for t in range(400):
        nb = int(rng.choice(cnt[use])); idx = rng.choice(len(r), nb, replace=False); rs.append(r[idx].sum() ** 2 / nb)
    sd_random = math.sqrt(np.mean(rs)); rat[s] = sd_blocks / sd_random
pr = {s: R['blocks']['ratio_blocks_random']['%s_k1' % s]['ratio'] for s in ('even', 'odd')}
check('consecutive blocks', all(abs(rat[s] / pr[s] - 1) < 0.25 for s in rat), 'k = 1 block/random scatter: own %.3f / %.3f, published %.3f / %.3f (even / odd)' % (rat['even'], rat['odd'], pr['even'], pr['odd']))
# ---------------------------------------------------------------- 5. small values
def beta(v, e0=0.03): v = v[v < e0]; b = len(v) / np.sum(np.log(e0 / v)); return b, b / math.sqrt(len(v))
b_all = beta(L); b_e = beta(L[D > 0]); b_o = beta(L[D < 0])
check('small values', abs(b_all[0] - R['small_values']['exponent']['vse']) < 1e-9 and abs(b_all[0] - 1.5) < 3 * b_all[1], 'exponent %.3f +- %.3f (even %.3f, odd %.3f), USp 3/2' % (b_all[0], b_all[1], b_e[0], b_o[0]))
small = sorted(zip(L.tolist(), D.tolist()))[:3]; ps = R['nonvanishing']['smallest'][:3]
check('smallest values', all(d == p[0] and abs(v - p[1]) < 1e-18 for (v, d), p in zip(small, ps)), ', '.join('D = %d: %.4g' % (d, v) for v, d in small))
# ---------------------------------------------------------------- 6. values at 10^10: own functional equation against PARI
De, Le = load('E10_E100'); Dp, Lp = load('E10_pari'); own = dict(zip(De.tolist(), Le.tolist())); rel = np.array([abs(own[d] - v) / v for d, v in zip(Dp.tolist(), Lp.tolist())])
check('values at 1e10', len(rel) == 12153 and rel.max() < 1e-9 and abs(rel.max() - R['checks']['e10_pari_vs_own']['max_rel']) < 1e-15, 'E10: %d values, largest relative difference %.2e, median %.1e' % (len(rel), rel.max(), np.median(rel)))
# 10^9: the family S9L = every fundamental discriminant with 10^9 < |D| <= 10^9 + 8*10^5 (own segmented sieve), its values against the 2,000 PARI values of the check sample,
# and the moments against the CFKRS polynomials (descriptive: no rule was recorded for them)
def fund_segment(lo, hi):
    def squarefree(a, b):
        nn = np.arange(a + 1, b + 1, dtype=np.int64); s_ok = np.ones(len(nn), bool)
        for q in range(2, int(b ** 0.5) + 1): s_ok[((-(a + 1)) % (q * q))::q * q] = False
        return nn, s_ok
    nn, s1 = squarefree(lo, hi); out = [nn[s1 & (nn % 4 == 1)], -nn[s1 & (nn % 4 == 3)]]
    mm, s2 = squarefree(lo // 4, hi // 4); s2 &= (4 * mm > lo); out += [4 * mm[s2 & ((mm % 4 == 2) | (mm % 4 == 3))], -4 * mm[s2 & ((mm % 4 == 1) | (mm % 4 == 2))]]
    return np.concatenate(out)
D9, L9 = load('S9L'); F9 = fund_segment(10 ** 9, 10 ** 9 + 8 * 10 ** 5); Dk, Lk = load('S9K_pari'); own9 = dict(zip(D9.tolist(), L9.tolist())); rel9 = np.array([abs(own9[d] - v) / v for d, v in zip(Dk.tolist(), Lk.tolist())])
S9 = R['weighted_density']['corrected']['S9L']; K9 = S9['central_values']['check_against_PARI']; d9 = 0.0
for s, sg in (('even', 1), ('odd', -1)):
    m = np.sign(D9) == sg; x = np.log(np.abs(D9[m]).astype(float))
    for k in range(1, 6): d9 = max(d9, abs((L9[m] ** k).sum() / model(k, s, x).sum() - 1 - S9['moments_vs_cfkrs']['%s_k%d' % (s, k)]['rel']))
check('values at 1e9', len(F9) == len(D9) == 486344 and np.array_equal(np.sort(F9), np.sort(D9)) and L9.min() > 0 and len(rel9) == 2000 and rel9.max() < 1e-9 and abs(rel9.max() - K9['max_relative_difference']) < 1e-15 and d9 < 1e-9,
      'S9L: own segmented sieve %d fundamental discriminants = published values; 2,000 PARI values: largest relative difference %.2e, median %.1e; moments k = 1..5 against CFKRS = published to %.1e (descriptive, largest %.2f errors)' % (len(F9), rel9.max(), np.median(rel9), d9, S9['moments_vs_cfkrs']['max_abs_z']))
# ---------------------------------------------------------------- 7. corrected weighted density (Section 4.6) from the published bins
WC = R['weighted_density']['corrected']; G = WC['groups']; u = np.array(R['weighted_density']['amplitude_shape_u']); dc = da = 0.0; pmin = 1.0
for k, g in G.items():
    d, p_, p0, e = (np.array(g[x]) for x in ('data', 'prediction', 'prediction_erroneous', 'error')); c2 = float(np.sum(((d - p_) / e) ** 2)); pmin = min(pmin, stats.chi2.sf(c2, 5))
    a = np.sum(u * (d - p_) / e ** 2) / np.sum(u * u / e ** 2); a0 = np.sum(u * (d - p0) / e ** 2) / np.sum(u * u / e ** 2)
    dc = max(dc, abs(c2 - g['chi2']) / max(1, g['chi2'])); da = max(da, abs(a - g['amplitude']), abs(a0 - g['amplitude_erroneous']))
amax = max(abs(g['amplitude']) / g['amplitude_error'] for k, g in G.items() if k.endswith('_even'))
n7 = sum(1 for d in fund if 7500000 < abs(d) <= 10000000); n7e = sum(1 for d in fund if 7500000 < d <= 10000000); S7 = WC['S7']; n9 = len(F9); n9e = int((F9 > 0).sum())
check('weighted density', dc < 1e-9 and da < 1e-9 and pmin > 1e-3 and amax < 3 and len(G) == 20 and n7 == S7['discriminants'] == S7['target_all_fundamental_7p5e6_1e7'] and n7e == G['S7_even']['n'] and n7 - n7e == G['S7_odd']['n'] and S7['missing'] == 0
      and n9 == S9['discriminants'] == S9['target_all_fundamental_1e9_1e9_plus_8e5'] and n9e == G['S9L_even']['n'] and n9 - n9e == G['S9L_odd']['n'] and S9['missing'] == 0 and S9['without_zeros'] == 0,
      '%d groups: chi2 and amplitudes = published to %.1e; smallest p %.2f; largest |even amplitude| %.1f errors (erroneous first implementation: C100 %.2f, S7 %.2f, S9L %.2f); S7 = all %d and S9L = all %d fundamental discriminants of their ranges (own sieves)' % (len(G), max(dc, da), pmin, amax, G['C100_even']['amplitude_erroneous'], G['S7_even']['amplitude_erroneous'], G['S9L_even']['amplitude_erroneous'], n7, n9))
EF = {F + '_' + zn: E for F, key in (('C100', 'explicit_formula_C100'), ('S7', 'explicit_formula_S7')) for zn, E in WC[key].items()}; ds = max(abs(E['delta_prime_side'] - E['parts']['gamma'] - sum(v for k, v in E['parts'].items() if k != 'gamma')) for E in EF.values())
zr = [E['zero_minus_prime_rule_tail'] / E['error_zero_minus_prime'] for E in EF.values()]
check('explicit formula', ds < 1e-12 and all(abs(z) < 3 for z in zr) and all(abs(E['delta_zero_side'] - E['delta_prime_side'] - E['zero_minus_prime_smooth_tail']) < 1e-12 for E in EF.values()),
      'C100 and S7: prime side = sum of its parts to %.0e; zero side - prime side %.1f / %.1f errors on C100 and %.1f / %.1f on S7 (even / odd)' % (ds, zr[0], zr[1], zr[2], zr[3]))
print('%d of %d checks passed' % (sum(ok), len(ok)))
sys.exit(0 if all(ok) else 1)
