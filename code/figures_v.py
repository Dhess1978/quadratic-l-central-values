#!/usr/bin/env python3
"""Figures of study V (light background, English and Czech), PNG 300 dpi and PDF, from data/results/study_v_results.json.
usage: python3 -I figures_v.py data/results/study_v_results.json figures/"""
import sys, json, math, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter
R = json.load(open(sys.argv[1])); OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'legend.fontsize': 8, 'axes.edgecolor': '#555555',
                     'axes.linewidth': 0.8, 'xtick.color': '#333333', 'ytick.color': '#333333', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': False})
INK, MUT, RAT, USP, EV, OD, KS, NUL = '#222222', '#8a8a8a', '#1f5fa8', '#d07a1a', '#1f5fa8', '#c0392b', '#7b5ea7', '#b0b0b0'
KCOL = ['#1f5fa8', '#2e8b57', '#d07a1a', '#c0392b', '#7b5ea7']
def cz(s): return s.replace('.', ',')
CZF = FuncFormatter(lambda v, p: ('%g' % round(v, 6)).replace('.', ',').replace('-', '−'))
def comma(fig):
    for ax in fig.axes:
        if ax.get_yscale() == 'linear' and isinstance(ax.yaxis.get_major_formatter(), ScalarFormatter): ax.yaxis.set_major_formatter(CZF)
        if ax.get_xscale() == 'linear' and isinstance(ax.xaxis.get_major_formatter(), ScalarFormatter): ax.xaxis.set_major_formatter(CZF)
def save(fig, name):
    if name.endswith('_cz'): comma(fig)
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300, bbox_inches='tight'); fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', metadata={'CreationDate': None}); plt.close(fig)
M = R['moments']; BL = R['blocks']; LV = R['log_values']; SV = R['small_values']; FZ0 = R['first_zero']; WD = R['weighted_density']
for en in (False, True):
    sfx = '_en' if en else '_cz'; w = lambda c, e: e if en else c; f = lambda x, d: ('%.*f' % (d, x)) if en else cz('%.*f' % (d, x))
    EVN, ODN = w('sudé (D > 0)', 'even (D > 0)'), w('liché (D < 0)', 'odd (D < 0)')
    # 1) moments k = 1..5 against CFKRS
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.7), gridspec_kw={'width_ratios': [1.35, 1]})
    ax = axs[0]; ax.axhspan(-1, 1, color=RAT, alpha=0.06, lw=0); ax.axhline(0, color=MUT, lw=0.8)
    for g in M['bands']:
        k = g['k']; xm = 0.5*(math.log10(max(g['lo'], 3)) + math.log10(g['hi'])); dx = (k - 3)*0.045 + (0.012 if g['s'] == 'odd' else -0.012)
        ax.plot(xm + dx, g['rel']/g['e'], 'o' if g['s'] == 'even' else 's', ms=3.6, color=KCOL[k - 1], mfc=KCOL[k - 1] if g['s'] == 'even' else 'white', mew=0.9)
    for key, g in M['at_1e8'].items():
        s, kk = key.split('_k'); k = int(kk); dx = (k - 3)*0.045 + (0.012 if s == 'odd' else -0.012)
        ax.plot(8.0 + dx, g['z'], 'o' if s == 'even' else 's', ms=3.6, color=KCOL[k - 1], mfc=KCOL[k - 1] if s == 'even' else 'white', mew=0.9)
    for k in range(1, 6): ax.plot([], [], 'o', color=KCOL[k - 1], ms=4, label='k = %d' % k)
    ax.plot([], [], 'o', color=INK, ms=4, label=EVN); ax.plot([], [], 's', color=INK, mfc='white', ms=4, label=ODN)
    ax.set_ylim(-2.2, 2.2); ax.set_xlabel(w('log₁₀|D| (střed pásma)', 'log₁₀|D| (centre of the band)')); ax.set_ylabel(w('(data − CFKRS) / chyba', '(data − CFKRS) / error'))
    ax.set_xticks([2, 3, 4, 5, 6, 7, 8]); ax.set_xlim(1.9, 8.4); ax.text(8.0, -1.95, w('D100b+c', 'D100b+c'), fontsize=7, ha='center')
    ax.legend(loc='lower left', ncol=4, fontsize=7); ax.set_title(w('Momenty k = 1–5 proti CFKRS: 7 pásem do 10⁷ a rodiny u 10⁸', 'Moments k = 1–5 against CFKRS: 7 bands up to 10⁷ and families at 10⁸'), fontsize=9)
    ax = axs[1]; ks = np.arange(1, 6)
    for s, mk, col in (('even', 'o', EV), ('odd', 's', OD)):
        full = [abs(M['total']['%s_k%d' % (s, k)]['rel']) for k in ks]; err = [M['total']['%s_k%d' % (s, k)]['err_emp_rel'] for k in ks]
        lead = [abs(M['leading_term_only_rel']['%s_k%d' % (s, k)]) for k in ks]
        ax.plot(ks, lead, mk + '--', color=KS, mfc=KS if s == 'even' else 'white', ms=4.5, lw=1, label=w('jen vedoucí člen (Keating–Snaith), ', 'leading term only (Keating–Snaith), ') + (EVN if s == 'even' else ODN))
        ax.plot(ks, full, mk + '-', color=col, mfc=col if s == 'even' else 'white', ms=4.5, lw=1, label=w('plný polynom CFKRS, ', 'full CFKRS polynomial, ') + (EVN if s == 'even' else ODN))
        ax.plot(ks, err, ':', color=col, lw=1)
    ax.plot([], [], ':', color=INK, lw=1, label=w('chyba z bloků', 'error from blocks'))
    ax.set_yscale('log'); ax.set_xticks(ks); ax.set_xlabel('k'); ax.set_ylabel(w('|relativní rozdíl od dat|', '|relative difference from the data|')); ax.legend(fontsize=6.6, loc='upper left')
    ax.set_title(w('Všechna D s |D| ≤ 10⁷ (6 079 285 L)', 'All D with |D| ≤ 10⁷ (6,079,285 L)'), fontsize=9); ax.set_ylim(1e-7, 1e6)
    save(fig, '01_moments' + sfx)
    # 2) consecutive blocks
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.6))
    ax = axs[0]
    for key, col, mk, ls in (('V_even_k1', EV, 'o', '-'), ('V_odd_k1', OD, 's', '-'), ('V_even_k2', EV, 'o', '--'), ('V_odd_k2', OD, 's', '--')):
        d = np.array(BL['width'][key]); ax.plot(d[:, 0], d[:, 1], mk + ls, color=col, ms=3.5, lw=1, mfc=col if 'even' in key else 'white', label=('|D| ≤ 10⁷, ' + (EVN if 'even' in key else ODN) + ', k = ' + key[-1]))
    for key, col, mk in (('D100_even_k1', EV, '^'), ('D100_odd_k1', OD, 'v')):
        d = np.array(BL['width'][key]); ax.plot(d[:, 0], d[:, 1], mk + ':', color=col, ms=4, lw=1, mfc='white', label='D100 (10⁸), ' + (EVN if 'even' in key else ODN) + ', k = 1')
    ax.axhline(1, color=MUT, lw=0.8, ls='--'); ax.set_xscale('log'); ax.set_ylim(0, 1.05); ax.set_xlabel(w('šířka bloku w (počet po sobě jdoucích |D|)', 'width of the block w (consecutive |D|)'))
    ax.set_ylabel(w('rozptyl součtů bloků / náhodné výběry', 'scatter of block sums / random samples')); ax.legend(fontsize=6.4, loc='upper right')
    ax.set_title(w('Souvislé bloky se vyruší (náhodné výběry = 1)', 'Consecutive blocks cancel (random samples = 1)'), fontsize=9)
    ax = axs[1]
    for s, col, mk in (('even', EV, 'o'), ('odd', OD, 's')):
        for k, ls in ((1, '-'), (2, '--')):
            ws = [1000, 3000, 10000, 30000]; al = [BL['scaling_alpha']['%s_k%d_w%d' % (s, k, x)]['alfa'] for x in ws]
            ax.plot(ws, al, mk + ls, color=col, ms=4, lw=1, mfc=col if s == 'even' else 'white', label=(EVN if s == 'even' else ODN) + ', k = %d' % k)
    ax.axhline(0.5, color=INK, lw=0.8, ls=':'); ax.axhline(0, color=INK, lw=0.8, ls=':')
    ax.text(1150, 0.47, w('w/√|D| (α = ½)', 'w/√|D| (α = ½)'), fontsize=7.5, va='top'); ax.text(1150, 0.02, w('jen w (α = 0)', 'w alone (α = 0)'), fontsize=7.5)
    ax.set_xscale('log'); ax.set_ylim(-0.05, 0.6); ax.set_xlabel(w('šířka bloku w u |D| ≈ 10⁸', 'width of the block w at |D| ≈ 10⁸')); ax.set_ylabel(w('efektivní exponent α ve w/|D|^α', 'effective exponent α in w/|D|^α'))
    ax.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=4); ax.set_title(w('Pásmo 10⁶–3·10⁶ proti 10⁸ + 10⁵ … 10⁸ + 2·10⁶', 'Band 10⁶–3·10⁶ against 10⁸ + 10⁵ … 10⁸ + 2·10⁶'), fontsize=9)
    save(fig, '02_blocks' + sfx)
    # 3) ln L against the hybrid model; fractional moments
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.6))
    ax = axs[0]
    for s, col, mk in (('even', EV, 'o'), ('odd', OD, 's')):
        G = [g for g in LV['bands'] if g['s'] == s]; x = [math.log10(math.sqrt(g['lo']*g['hi'])) for g in G]
        ax.plot(x, [g['data'][1]/g['pred'][1] - 1 for g in G], mk + '-', color=col, ms=4, lw=1, mfc=col if s == 'even' else 'white', label=w('rozptyl, ', 'variance, ') + (EVN if s == 'even' else ODN))
        ax.plot(x, [g['data'][2]/g['pred'][2] - 1 for g in G], mk + '--', color=col, ms=4, lw=1, mfc='white', label=w('třetí kumulant, ', 'third cumulant, ') + (EVN if s == 'even' else ODN))
        ax.plot(x, [g['data'][0] - g['pred'][0] for g in G], mk + ':', color=col, ms=3, lw=1, mfc='white', alpha=0.7, label=w('střední hodnota (absolutně), ', 'mean (absolute), ') + (EVN if s == 'even' else ODN))
    ax.axhline(0, color=MUT, lw=0.8); ax.set_xlabel('log₁₀|D|'); ax.set_ylabel(w('data / hybrid − 1', 'data / hybrid − 1')); ax.legend(fontsize=6.3, loc='lower right', ncol=1)
    ax.set_title(w('ln L(½) proti hybridu (N z prvního momentu)', 'ln L(½) against the hybrid (N from the first moment)'), fontsize=9)
    ax = axs[1]; P = LV['fractional_trend']['pasma']; FT = LV['fractional_trend']['fit']
    for s, col in (('even', EV), ('odd', OD)):
        for sx, mk in (('0.5', 'o'), ('-0.5', 'v')):
            G = [p for p in P if p['s'] == s]; xx = np.array([1/math.log(math.sqrt(p['lo']*p['hi'])) for p in G])
            ax.errorbar(xx, [100*p['rel'][sx] for p in G], yerr=[200*p['e'][sx] for p in G], fmt=mk, color=col, ms=4, mfc=col if sx == '0.5' else 'white', capsize=2, lw=0.9)
            ft = FT['%s_%s' % (s, sx)]; xs = np.linspace(0, 0.09, 50); ax.plot(xs, 100*(ft['a'] + ft['b']*xs), '-', color=col, lw=0.8, alpha=0.6)
    ax.plot([], [], 'o', color=INK, ms=4, label='E[L^½]'); ax.plot([], [], 'v', color=INK, mfc='white', ms=4, label='E[L^−½]'); ax.plot([], [], '-', color=EV, label=EVN); ax.plot([], [], '-', color=OD, label=ODN)
    ax.axhline(0, color=MUT, lw=0.8); ax.set_xlim(-0.002, 0.092); ax.set_xlabel('1 / ln|D|'); ax.set_ylabel(w('data / hybrid − 1 (%), ±2 chyby', 'data / hybrid − 1 (%), ±2 errors')); ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.set_title(w('Zlomkové momenty: odchylka klesá jako 1/ln|D|', 'Fractional moments: the deviation falls like 1/ln|D|'), fontsize=9)
    save(fig, '03_log_values' + sfx)
    # 4) small values
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.6))
    ax = axs[0]; C = np.array(SV['counts']); eps = C[:, 0]
    for j, s, col, mk in ((1, 'even', EV, 'o'), (2, 'odd', OD, 's')): ax.plot(eps, C[:, j], mk + '-', color=col, ms=4, lw=1, mfc=col if s == 'even' else 'white', label=EVN if s == 'even' else ODN)
    e0 = 0.03; ref = 0.5*(C[3, 1] + C[3, 2]); xs = np.logspace(-4, 0, 30); ax.plot(xs, ref*(xs/e0)**1.5, ':', color=INK, lw=1, label=w('sklon 3/2 (USp)', 'slope 3/2 (USp)'))
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('ε'); ax.set_ylabel(w('počet D s L(½) < ε (|D| ≤ 10⁷)', 'number of D with L(½) < ε (|D| ≤ 10⁷)')); ax.legend(fontsize=7.5, loc='upper left')
    ax.set_title(w('Malé hodnoty: exponent %s ± %s (USp 3/2)' % (f(SV['exponent']['vse'], 3), f(0.015, 3)), 'Small values: exponent %s ± %s (USp 3/2)' % (f(SV['exponent']['vse'], 3), f(0.015, 3))), fontsize=9)
    ax = axs[1]
    for s, col, mk in (('even', EV, 'o'), ('odd', OD, 's')):
        G = [g for g in SV['bands'] if g['s'] == s and g['lo'] >= 1e4]; x = np.array([math.log10(math.sqrt(g['lo']*g['hi'])) for g in G])
        for j, (ls, eps_) in enumerate((('-', '0.1'), ('--', '0.03'))):
            r = [g['obs'][j]/g['pred'][j] for g in G]; e = [math.sqrt(max(g['obs'][j], 1))/g['pred'][j] for g in G]
            ax.errorbar(x + (0.03 if s == 'odd' else -0.03) + 0.015*j, r, yerr=e, fmt=mk + ls, color=col, ms=4, lw=0.9, mfc=col if s == 'even' else 'white', capsize=2, alpha=1 - 0.35*j,
                        label=(EVN if s == 'even' else ODN) + ', ε = ' + (eps_ if en else cz(eps_)))
    ax.axhline(1, color=MUT, lw=0.8); ax.set_ylim(0.6, 1.25); ax.set_xlabel('log₁₀|D|'); ax.set_ylabel(w('pozorováno / hybrid', 'observed / hybrid')); ax.legend(fontsize=7, loc='upper right', ncol=2)
    ax.set_title(w('Počet L(½) < ε proti hybridu (±1 chyba)', 'Number of L(½) < ε against the hybrid (±1 error)'), fontsize=9)
    save(fig, '04_small_values' + sfx)
    # 5) central value and first zero
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.6)); FAM = ['B10', 'C100', 'R10', 'D100']; X = np.arange(4)
    ax = axs[0]; CT = FZ0['corr_theta']
    for s, col, mk, dx in (('even', EV, 'o', -0.12), ('odd', OD, 's', 0.12)):
        ax.plot(X + dx, [CT['%s_%s' % (F, s)]['corr'] for F in FAM], mk, color=col, ms=5, mfc=col if s == 'even' else 'white', label=w('data, ', 'data, ') + (EVN if s == 'even' else ODN))
        ax.plot(X + dx, [CT['%s_%s' % (F, s)]['corr_usp'] for F in FAM], '_', color=USP, ms=12, mew=2)
    ax.plot([], [], '_', color=USP, ms=12, mew=2, label=w('USp(2N), N z prvního momentu', 'USp(2N), N from the first moment'))
    ax.set_xticks(X); ax.set_xticklabels(FAM); ax.set_ylabel('corr(ln L(½), x₁)'); ax.legend(fontsize=7.5, loc='upper right'); ax.set_ylim(0.78, 0.92)
    ax.set_title(w('Hodnota ve středu a první nula: vazba silnější než USp', 'Central value and first zero: the coupling is stronger than USp'), fontsize=9)
    ax = axs[1]; NE = FZ0['n_eff']
    for s, col, mk, dx in (('even', EV, 'o', -0.12), ('odd', OD, 's', 0.12)):
        ax.plot(X + dx, [NE['%s_%s' % (F, s)]['d1']/NE['%s_%s' % (F, s)]['d1_usp'] for F in FAM], mk, color=col, ms=5, mfc=col if s == 'even' else 'white', label=w('posun při vážení L(½) / USp(2N_ef), ', 'shift under weighting by L(½) / USp(2N_eff), ') + (EVN if s == 'even' else ODN))
        ax.plot(X + dx, [NE['%s_%s' % (F, s)]['mx1']/NE['%s_%s' % (F, s)]['mx1_usp'] for F in FAM], mk, color=col, ms=4, mfc='white', alpha=0.6, label=w('střední x₁ / USp(2N_ef), ', 'mean x₁ / USp(2N_eff), ') + (EVN if s == 'even' else ODN))
    ax.axhline(1, color=MUT, lw=0.8); ax.set_xticks(X); ax.set_xticklabels(FAM); ax.set_ylim(0.97, 1.18); ax.set_ylabel(w('data / USp(2N_ef)', 'data / USp(2N_eff)')); ax.legend(fontsize=6.6, loc='center')
    ax.set_title(w('N_ef z domněnky poměrů: střední x₁ sedí, posun ne', 'N_eff from the ratios conjecture: mean x₁ agrees, the shift does not'), fontsize=9)
    save(fig, '05_first_zero' + sfx)
    # 6) weighted one-level density (corrected prediction, lab 34.85)
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.7)); xc = np.array([0.1, 0.3, 0.5, 0.7, 0.9]); C = WD['corrected']['groups']; G = WD['groups']
    ax = axs[0]
    for s, col, mk, dx in (('even', EV, 'o', -0.012), ('odd', OD, 's', 0.012)):
        g = C['C100_' + s]; ax.plot(xc, g['prediction'], '-', color=col, lw=1.2); ax.plot(xc, G['C100_' + s]['usp'], '--', color=USP, lw=0.9)
        ax.errorbar(xc + dx, g['data'], yerr=2*np.array(g['error']), fmt=mk, color=col, ms=4.5, mfc=col if s == 'even' else 'white', capsize=2, lw=0.9, label=(EVN if s == 'even' else ODN))
    ax.plot([], [], '-', color=INK, lw=1.2, label=w('domněnka poměrů s vahou L(½) (Fazzari)', 'ratios conjecture weighted by L(½) (Fazzari)')); ax.plot([], [], '--', color=USP, lw=0.9, label='USp(2N_x₁)')
    ax.axhline(0, color=MUT, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel(w('rozvinutá výška x (koše 0,2)', 'unfolded height x (bins of 0.2)')); ax.set_ylabel(w('vážená − nevážená hustota', 'weighted − unweighted density'))
    ax.legend(fontsize=7, loc='upper left'); ax.set_title(w('C100 (10⁶–2·10⁶, 607 941 L)', 'C100 (10⁶–2·10⁶, 607,941 L)'), fontsize=9)
    ax = axs[1]; Y6 = (float(os.environ.get('FIG6_YMIN', '-0.088')), float(os.environ.get('FIG6_YMAX', '0.042')))
    SETS = [('C100 (10⁶)', C['C100_even'], C['C100_odd'], -0.056), ('S7 (10⁷)', C['S7_even'], C['S7_odd'], -0.028), ('10⁸', C['D100_D100b_D100c_even'], C['D100_D100b_D100c_odd'], 0.0), ('S9L (10⁹)', C['S9L_even'], C['S9L_odd'], 0.028), ('E (10¹⁰)', C['E_even'], C['E_odd'], 0.056)]
    for (nm, ge, go, dx), col in zip(SETS, (EV, '#2a9d8f', INK, KS, USP)):
        for g, mk, ddx in ((ge, 'o', -0.006), (go, 's', 0.006)):
            r = np.array(g['data']) - np.array(g['prediction']); e = np.array(g['error'])
            ax.errorbar(xc + dx + ddx, r, yerr=2*e, fmt=mk + '-', color=col, ms=4, lw=0.9, mfc=col if mk == 'o' else 'white', capsize=2, alpha=1 if mk == 'o' else 0.75,
                        label=nm + ', ' + (w('sudé', 'even') if mk == 'o' else w('liché', 'odd')) + ' (χ² %s)' % f(g['chi2'], 1))
    ax.axhline(0, color=MUT, lw=0.8); ax.set_xlim(0, 1); ax.set_ylim(Y6[0], Y6[1]); ax.set_xlabel(w('rozvinutá výška x (koše 0,2)', 'unfolded height x (bins of 0.2)')); ax.set_ylabel(w('data − předpověď, ±2 chyby', 'data − prediction, ±2 errors'))
    ax.legend(fontsize=6.3, loc='lower left', ncol=2); ax.set_title(w('Zbytek proti domněnce poměrů: obě znaménka sedí', 'Residual against the ratios conjecture: both signs agree'), fontsize=9)
    save(fig, '06_weighted_density' + sfx)
    # 7) the apparent even term of the first implementation and the corrected prediction
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.8)); FF = ['B10', 'C100', 'R10', 'S7', 'D100c', 'D100_D100b_D100c', 'S9L', 'E']; LAB = {'D100_D100b_D100c': w('10⁸', '10⁸'), 'D100c': 'D100c', 'E': 'E', 'S9L': w('S9L (10⁹)', 'S9L (10⁹)')}; OFS = {'R10': (-19, 3), 'S7': (5, 5)}
    ax = axs[0]
    for F in FF:
        q, o = C[F + '_even'], C[F + '_odd']; dx = -0.45 if F == 'D100c' else 0
        ax.plot(q['l'] + dx - 0.1, q['amplitude_erroneous'], 'o', color=MUT, ms=5, mfc='white')
        ax.errorbar(q['l'] + dx, q['amplitude'], yerr=2*q['amplitude_error'], fmt='o', color=EV, ms=5, capsize=2.5, lw=1)
        ax.errorbar(q['l'] + dx + 0.12, o['amplitude'], yerr=2*o['amplitude_error'], fmt='s', color=OD, ms=3.5, mfc='white', capsize=2, lw=0.8, alpha=0.8)
        if F == 'D100c': ax.annotate('D100c', (q['l'] + dx, q['amplitude'] - 2*q['amplitude_error']), textcoords='offset points', xytext=(-14, -12), fontsize=7)
        else: ax.annotate(LAB.get(F, F), (q['l'] + dx, max(q['amplitude_erroneous'], q['amplitude'] + 2*q['amplitude_error'])), textcoords='offset points', xytext=OFS.get(F, (3, 5)), fontsize=7)
    ax.plot([], [], 'o', color=MUT, mfc='white', ms=5, label=w('sudé D, chybná první implementace', 'even D, erroneous first implementation'))
    ax.errorbar([], [], yerr=[], fmt='o', color=EV, ms=5, label=w('sudé D, opravená předpověď (±2 chyby)', 'even D, corrected prediction (±2 errors)'))
    ax.errorbar([], [], yerr=[], fmt='s', color=OD, ms=3.5, mfc='white', label=w('liché D, opravená předpověď', 'odd D, corrected prediction'))
    ax.axhline(0, color=MUT, lw=0.8); ax.set_xlim(8.5, 23.5); ax.set_ylim(-0.7, 2.7); ax.set_xlabel('l = ln(|D|/π)'); ax.set_ylabel(w('amplituda na tvaru dřívější odchylky C100', 'amplitude on the shape of the earlier C100 deviation'))
    ax.legend(fontsize=6.8, loc='upper right'); ax.set_title(w('Zdánlivý člen sudých D: chyba první implementace', 'Apparent term of even D: error of the first implementation'), fontsize=8.5)
    ax = axs[1]; FA = ['B10', 'C100', 'R10', 'S7', 'D100', 'D100_D100b', 'D100c', 'D100_D100b_D100c', 'S9L', 'E']; X = np.arange(len(FA))
    for s, col, mk, dx in (('even', EV, 'o', -0.15), ('odd', OD, 's', 0.15)):
        ax.plot(X + dx, [C[F + '_' + s]['chi2_erroneous'] for F in FA], mk, color=col, ms=5, mfc='white', alpha=0.7)
        ax.plot(X + dx, [C[F + '_' + s]['chi2'] for F in FA], mk, color=col, ms=5, label=(EVN if s == 'even' else ODN) + w(', opravená', ', corrected'))
    ax.plot([], [], 'o', color=MUT, mfc='white', ms=5, label=w('prázdné: chybná první implementace', 'open: erroneous first implementation')); ax.axhline(20.52, color=USP, lw=0.9, ls='--', label=w('p = 0,001', 'p = 0.001'))
    ax.set_yscale('log'); ax.set_ylim(0.1, 1000); ax.set_xticks(X); ax.set_xticklabels([{'D100_D100b': 'D100+b', 'D100_D100b_D100c': '10⁸', 'D100c': 'D100c'}.get(F, F) for F in FA], fontsize=7, rotation=25, ha='right', rotation_mode='anchor')
    ax.set_ylabel(w('χ² na 5 koších', 'χ² for 5 bins')); ax.legend(fontsize=6.8, loc='upper right'); ax.set_title(w('χ² vážené hustoty, chybná a opravená předpověď', 'χ² of the weighted density, erroneous and corrected'), fontsize=8.5)
    save(fig, '07_correction' + sfx)
    # 8) explicit formula: prime side by parts and zero side, C100 and S7 (Section 4.6)
    fig, axs = plt.subplots(1, 3, figsize=(10.6, 3.8), gridspec_kw={'width_ratios': [1, 1, 0.78]}); XF = {'C100': WD['corrected']['explicit_formula_C100'], 'S7': WD['corrected']['explicit_formula_S7']}
    PK = ['p2', 'odd_p_le_sqrt', 'odd_p_gt_sqrt', 'prime_powers']; PL = ['p = 2', w('lichá p\n≤ √(|D|/π)', 'odd p\n≤ √(|D|/π)'), w('lichá p\n> √(|D|/π)', 'odd p\n> √(|D|/π)'), w('mocniny\nprvočísel', 'prime\npowers'), w('strana prvočísel\n(součet)', 'prime side\n(sum)')]
    for ax, (F, nL) in zip(axs[:2], (('C100', w('607 941 L', '607,941 L')), ('S7', w('1 519 830 L', '1,519,830 L')))):
        for s_, col, dx in (('even', EV, -0.2), ('odd', OD, 0.2)):
            E = XF[F][s_]; v = [E['parts'][k] for k in PK] + [E['delta_prime_side']]; xs = np.arange(len(v)) + dx
            ax.bar(xs, v, width=0.36, color=col, alpha=0.82, lw=0, label=(EVN if s_ == 'even' else ODN))
            ax.plot([4 + dx], [E['delta_zero_side']], 'D', color=INK, ms=4.5, mfc='white', mew=1.1, label=w('strana nul', 'zero side') if s_ == 'even' else None, zorder=5)
            ax.annotate(f(v[4], 3).replace('-', '−'), (4 + dx + (-0.22 if s_ == 'even' else 0.22), v[4]), ha='right' if s_ == 'even' else 'left', va='center', fontsize=6.8, color=INK)
        ax.axhline(0, color=MUT, lw=0.8); ax.set_xticks(np.arange(5)); ax.set_xticklabels(PL, fontsize=6.8); ax.set_ylim(-0.56, 0.09); ax.set_xlim(-0.6, 5.35); ax.grid(axis='x', visible=False)
        ax.set_title(F + ' (' + nL + ')', fontsize=9); ax.legend(fontsize=7, loc='lower left')
        if F == 'C100': ax.set_ylabel(w('vážené − nevážené (příspěvek k Δ)', 'weighted − unweighted (contribution to Δ)'))
    ax = axs[2]; ax.axhspan(-3, 3, color=RAT, alpha=0.08, lw=0); ax.axhline(0, color=MUT, lw=0.8); k_ = 0; tl = []
    for F in ('C100', 'S7'):
        for s_, col in (('even', EV), ('odd', OD)):
            E = XF[F][s_]; z0 = E['zero_minus_prime_smooth_tail']/E['error_zero_minus_prime']; z1 = E['z_rule_tail']
            ax.plot([k_], [z0], 'o', color=col, ms=5.5, mfc='white', mew=1.1); ax.plot([k_], [z1], 'o', color=col, ms=5.5); ax.plot([k_, k_], [z0, z1], '-', color=col, lw=0.7, alpha=0.5)
            ax.annotate(f(z1, 1).replace('-', '−'), (k_, z1), textcoords='offset points', xytext=(7, -2), fontsize=6.8, color=INK); ax.annotate(f(z0, 1).replace('-', '−'), (k_, z0), textcoords='offset points', xytext=(7, -2), fontsize=6.8, color=MUT)
            tl.append(F + '\n' + (w('sudé', 'even') if s_ == 'even' else w('liché', 'odd'))); k_ += 1
    ax.plot([], [], 'o', color=INK, ms=5.5, label=w('ocas z hustot domněnky poměrů (pravidlo)', 'tail from the ratios-conjecture densities (rule)')); ax.plot([], [], 'o', color=INK, ms=5.5, mfc='white', label=w('ocas jen z hladké hustoty', 'tail from the smooth density alone'))
    ax.set_xticks(range(4)); ax.set_xticklabels(tl, fontsize=7); ax.set_xlim(-0.6, 3.9); ax.set_ylim(-80, 12); ax.set_ylabel(w('(strana nul − strana prvočísel) / párová chyba', '(zero side − prime side) / paired error')); ax.grid(axis='x', visible=False)
    ax.legend(fontsize=6.0, loc='lower center'); ax.set_title(w('Vyrovnání obou stran (pásmo ±3)', 'Balance of the two sides (band ±3)'), fontsize=9)
    fig.tight_layout(); save(fig, '08_explicit_formula' + sfx)
    # 9) appendix: Ford–Zaharescu
    FZ = R['ford_zaharescu']; B6 = ['Z6a', 'Z6b', 'Z6c', 'Z6d']; NB = FZ['bins']; Nn = sum(FZ['blocks'][b]['N'] for b in B6); wv = sum(FZ['blocks'][b]['N']*FZ['blocks'][b]['w'] for b in B6)/Nn
    def dens(u, a):
        p = {'p2': 2, 'p7': 7}.get(a); s = np.ones_like(u)
        if p:
            for k in range(1, FZ['KMAX'] + 1): s -= 2*wv*math.log(p)/p**(k/2)*np.cos(2*np.pi*k*u)
        return s
    fig, ax = plt.subplots(figsize=(7.4, 3.5)); us = np.linspace(0, 1, 400)
    for a, col, mk, nm in (('p2', RAT, 'o', 'α = ln 2/2π'), ('p7', USP, 'o', 'α = ln 7/2π'), ('odm2', MUT, 's', w('α = √2/2π (kontrola)', 'α = √2/2π (control)'))):
        O = np.zeros(NB); s2 = np.zeros(NB)
        for b in B6:
            x = FZ['blocks'][b]['alphas'][a]; O += np.array(x['O']); s2 += np.maximum(np.array(x['E']), np.array(x['sdn'])**2)
        d = O/(Nn/NB); e = np.sqrt(s2)/(Nn/NB); xc = (np.arange(NB) + 0.5)/NB
        ax.plot(us, dens(us, a), '-', color=col, lw=1.1); ax.errorbar(xc, d, yerr=2*e, fmt=mk, color=col, ms=2.4, mfc=col if a != 'odm2' else 'white', capsize=0, lw=0.6, label=nm)
    ax.set_ylim(0.62, 1.25); ax.set_xlim(0, 1); ax.set_xlabel('u = {α·γₙ}'); ax.set_ylabel(w('hustota (rovnoměrná = 1)', 'density (uniform = 1)')); ax.legend(fontsize=7.5, loc='lower center', ncol=3)
    ax.set_title(w('Ford–Zaharescu na 2 001 052 nulách ζ: čáry = Landauův vzorec bez volného parametru', 'Ford–Zaharescu on 2,001,052 zeros of ζ: lines = Landau formula without a free parameter'), fontsize=9)
    save(fig, '09_ford_zaharescu' + sfx)
print('figures written to', OUT)
