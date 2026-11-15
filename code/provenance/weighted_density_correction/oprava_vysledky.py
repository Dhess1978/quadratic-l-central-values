"""34.85: vážená hustota s opravenou předpovědí (nezávislá implementace) — χ², p, zbytková amplituda na tvaru dosavadní odchylky C100 sudých,
shoda opraveného kódu 34.79 s nezávislou implementací po rodinách."""
import json, math, numpy as np
from scipy import stats
import importlib.util, sys
R = json.load(open('/home/claude/pub5/repo/data/results/study_v_results.json'))['weighted_density']; N = json.load(open('/home/claude/br85/data/nez_porovnani.json'))
spec = importlib.util.spec_from_file_location('oq', '/home/claude/br85/py/oprava_q4.py')
src = open('/home/claude/br85/py/oprava_q4.py').read().split('PC0 = ')[0]; oq = {}; exec(src, oq)       # opravené kusy 34.79 (bez tisku)
s = np.array(R['amplitude_shape_u'])      # dosavadní tvar odchylky C100 sudých (data − chybná předpověď)
out = {'priciny': 'fazzari.py (34.79): ve členu s oběma výměnami α → −α, β → −β Γ-podíl gS(½+α+β) místo gS(½+α)·gS(½+β); táž chyba ve fazzari_tridy.py',
       'skupiny': {}}
for key in ['B10', 'C100', 'R10', 'D100', 'D100_D100b', 'D100c', 'D100_D100b_D100c', 'E']:
    for zn in ('even', 'odd'):
        g = R['groups'][key + '_' + zn] if key in ('B10', 'C100', 'R10', 'D100') else R[key][zn]
        d = np.array(g['rozdil']); eb = np.array(g['chyba'])
        if key in ('B10', 'C100', 'R10', 'D100'): e = np.maximum(eb, np.array(R['null'][key + '_' + zn]['nul_sd']))
        else: e = eb
        n = N[key + '_' + zn]; pr = np.array(n['pred_nez']); old = np.array(g['pred'])
        r = d - pr; c2 = float(np.sum((r / e) ** 2)); p = float(stats.chi2.sf(c2, 5)); c2b = float(np.sum((r / eb) ** 2))
        c2o = float(np.sum(((d - old) / e) ** 2)); po = float(stats.chi2.sf(c2o, 5))
        wgt = s / e ** 2; a = float(np.sum(r * wgt) / np.sum(s * wgt)); ea = float(1 / math.sqrt(np.sum(s * s / e ** 2)))
        ao = float(np.sum((d - old) * wgt) / np.sum(s * wgt))
        out['skupiny'][key + '_' + zn] = dict(n=n['n'], l=g.get('l'), rozdil=d.tolist(), chyba=e.tolist(), pred_opravena=pr.tolist(), pred_chybna=old.tolist(), chi2=c2, p=p, chi2_bootstrap=c2b,
                                              chi2_chybna=c2o, p_chybna=po, amplituda=a, amplituda_chyba=ea, amplituda_chybna=ao)
        print('%-22s χ² %6.2f p %.3f (bootstrap χ² %5.2f) | chybná χ² %7.1f | amplituda %+.3f ± %.3f (chybná %+.3f)' % (key + '_' + zn, c2, p, c2b, c2o, a, ea, ao))
# shoda opraveného kódu 34.79 s nezávislou implementací (stejné binování jako 34.79: 20 košů l)
CV = '/home/claude/pub5/repo/data/central_values/'
V = np.concatenate([np.load(CV + f)['D'] for f in sorted(__import__('os').listdir(CV)) if f.startswith('V_part')])
Dsets = {'C100': V[(np.abs(V) > 10 ** 6) & (np.abs(V) <= 2 * 10 ** 6)], 'E': np.load(CV + 'E10_E100.npz')['D'], 'B10': V[(np.abs(V) > 10 ** 4) & (np.abs(V) <= 10 ** 5)]}
EDG = np.linspace(0, 1, 6); BW = 0.2; shoda = {}
for F, Ds in Dsets.items():
    for zn, sg, a in (('even', 1, 0), ('odd', -1, 1)):
        lv = np.log(np.abs(Ds[np.sign(Ds) == sg]) / math.pi); ed = np.linspace(lv.min(), lv.max() + 1e-9, 21); cnt, _ = np.histogram(lv, ed)
        cen = [lv[(lv >= ed[i]) & (lv < ed[i + 1])].mean() for i in range(20) if cnt[i] > 0]; wts = cnt[cnt > 0].astype(float)
        num = np.zeros(5); den = 0.0
        for w, lc in zip(wts, cen):
            M, cw = oq['vahova'](lc, a, oq['PC']); num += w * (M.real * BW + oq['kose'](lc, a, cw, EDG)); den += w * M.real * BW
        vaz = num / den; vn = np.array(N[F + '_' + zn]['vazena_nez']); e = np.array(out['skupiny'][F + '_' + zn]['chyba'])
        shoda[F + '_' + zn] = float(np.max(np.abs(vaz - vn) / e)); print('shoda opravený 34.79 × nezávislá %-10s max %.4f chyby (abs %.1e)' % (F + '_' + zn, shoda[F + '_' + zn], np.max(np.abs(vaz - vn))))
out['shoda_opraveny_kod_s_nezavislym'] = shoda
json.dump(out, open('/home/claude/br85/data/oprava_vysledky.json', 'w'), indent=1)
