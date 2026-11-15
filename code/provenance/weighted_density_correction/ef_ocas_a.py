"""34.85a: ocas strany nul nad spočtenými nulami (x > 4,35) z hustot domněnky poměrů (vážené a nevážené, opravená předpověď) — podle pravidla odv-lv-ef;
zobecnění br85/py/ef_ocas.py na libovolnou rodinu. usage: python3 -I ef_ocas_a.py RODINA ADRESAR_DAT   (čte ef_RODINA_{even,odd}.json + _D.npz, píše ef_ocas_RODINA.json)"""
import sys, json, math, numpy as np, scipy.special as sp
sys.path.insert(0, '/home/claude/br85/py'); import nez_porovnani as NP
F, B = sys.argv[1], sys.argv[2].rstrip('/') + '/'
phi = lambda u: np.sinc(u / 2) ** 4; out = {}
for zn, a in (('even', 0), ('odd', 1)):
    z = np.load(B + 'ef_%s_%s_D.npz' % (F, zn)); l = np.log(np.abs(z['D']) / math.pi); L = z['L']
    ed = np.linspace(l.min(), l.max() + 1e-9, 41); idx = np.digitize(l, ed) - 1; vw = nm = vu = nw = 0.0; vw2 = vu2 = 0.0
    for i in range(40):
        m = idx == i
        if not m.any(): continue
        w = m.sum(); lc = l[m].mean(); Mv, r0, s_w, s_u = NP.hustoty(lc, a); TG = NP.TG
        th = 0.5 * TG * lc + np.imag(sp.loggamma((0.5 + a + 1j * TG) / 2)); oc = th / math.pi > 4.35; ph = phi(TG * lc / (2 * math.pi)) * oc
        Lm = L[m].mean()      # váha skupiny podle naměřeného ⟨L⟩ (jako v datech)
        vw += Lm * w * np.sum(ph * (r0 + s_w / Mv)) * NP.DT; nm += Lm * w; vu += w * np.sum(ph * (r0 + s_u)) * NP.DT; nw += w
        vw2 += Lm * w * np.sum(ph * r0) * NP.DT; vu2 += w * np.sum(ph * r0) * NP.DT
    d = 2 * (vw / nm - vu / nw); d0 = 2 * (vw2 / nm - vu2 / nw)
    J = json.load(open(B + 'ef_%s_%s.json' % (F, zn))); zp = J['delta']['Z-P']; e = J['chyba']['Z-P']; zp2 = zp - J['delta']['ocas'] + d
    out[zn] = dict(ocas_pravidlo=float(d), ocas_hladky_model=float(d0), ocas_v_ef=float(J['delta']['ocas']), Z_P_hladky=float(zp), Z_P_pravidlo=float(zp2), chyba=float(e), z_hladky=float(zp / e), z_pravidlo=float(zp2 / e))
    print('%s %s: ocas z domněnky poměrů Δ = %.3e (hladký model %.3e, v ef.py %.3e); Z − P %.2e ± %.1e (%.1f chyby) → s ocasem podle pravidla %.2e (%.1f chyby)' % (F, zn, d, d0, J['delta']['ocas'], zp, e, zp / e, zp2, zp2 / e))
json.dump(out, open(B + 'ef_ocas_%s.json' % F, 'w'), indent=1)
