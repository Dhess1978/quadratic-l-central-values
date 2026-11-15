"""ocas strany nul nad spočtenými nulami (x > 4,35) z hustot domněnky poměrů (vážené a nevážené, opravená předpověď) — podle pravidla odv-lv-ef"""
import sys, json, math, numpy as np, scipy.special as sp
sys.path.insert(0, '/home/claude/br85/py'); import nez_porovnani as NP
phi = lambda u: np.sinc(u / 2) ** 4
for zn, a in (('even', 0), ('odd', 1)):
    z = np.load('/home/claude/br85/data/ef2_C100_%s_D.npz' % zn); l = np.log(np.abs(z['D']) / math.pi); L = z['L']
    ed = np.linspace(l.min(), l.max() + 1e-9, 41); idx = np.digitize(l, ed) - 1; vw = nm = vu = nw = 0.0; vw2 = nm2 = vu2 = 0.0
    for i in range(40):
        m = idx == i
        if not m.any(): continue
        w = m.sum(); lc = l[m].mean(); Mv, r0, s_w, s_u = NP.hustoty(lc, a); TG = NP.TG
        th = 0.5 * TG * lc + np.imag(sp.loggamma((0.5 + a + 1j * TG) / 2)); oc = th / math.pi > 4.35; ph = phi(TG * lc / (2 * math.pi)) * oc
        Lm = L[m].mean()      # váha skupiny podle naměřeného ⟨L⟩ (jako v datech)
        vw += Lm * w * np.sum(ph * (r0 + s_w / Mv)) * NP.DT; nm += Lm * w; vu += w * np.sum(ph * (r0 + s_u)) * NP.DT; nw += w
        vw2 += Lm * w * np.sum(ph * r0) * NP.DT; vu2 += w * np.sum(ph * r0) * NP.DT
    d = 2 * (vw / nm - vu / nw); d0 = 2 * (vw2 / nm - vu2 / nw)
    J = json.load(open('/home/claude/br85/data/ef2_C100_%s.json' % zn)); zp = J['delta']['Z-P']; e = J['chyba']['Z-P']
    print('%s: ocas z domněnky poměrů Δ = %.2e (hladký model %.2e); Z − P původně %.2e ± %.1e → s ocasem podle pravidla %.2e (%.1f chyby)' % (zn, d, d0, zp, e, zp - d0 + d, (zp - d0 + d) / e))
