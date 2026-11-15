"""34.79 odv-lr-b34 podle PRAVIDLA_3479: člen x^{3/4} vyhlazeného třetího momentu v podrodině D = 8d (d liché bez čtverců, D > 0) proti Diaconu–Whitehead (Theorem 1.1–1.2, Res ≈ −0,0034).
R(x) = Σ_d (L³ − Q₃^{8d}(ln 8d))·W(d/x), b = R/(Ŵ(¾)·x^{3/4}); okna x = 1,25·10⁶/2^j (j = 0…5, nosiče [x/2, x] disjunktní); chyba z rozptylu oken."""
import json, math, numpy as np
from scipy.integrate import quad
V = [np.load('/home/claude/br76/data/%s.npz' % f) for f in ('V5', 'V10')]; D = np.concatenate([v['D'] for v in V]).astype(np.int64); L = np.concatenate([v['L'] for v in V])
m = (D > 0) & (D % 16 == 8); D8 = D[m]; L8 = L[m]; d = D8 // 8; assert np.all(d % 2 == 1)
c = json.load(open('/home/claude/br79/data/q3sub.json'))['sub']; Q = np.polyval(np.array(c)[::-1], np.log(D8.astype(float))); r = L8 ** 3 - Q
W = lambda u: np.where((u > 0.5) & (u < 1), (1 - (4 * u - 3) ** 2) ** 4, 0.0); Wh = quad(lambda u: float(W(np.array(u))) * u ** -0.25, 0.5, 1)[0]
RES = -0.0034; out = dict(n=int(len(d)), Wh=Wh, res=RES, okna=[])
for j in range(6):
    x = 1.25e6 / 2 ** j; w = W(d / x); R = float(np.sum(r * w)); M = float(np.sum(Q * w)); Dm = float(np.sum(L8 ** 3 * w)); eiid = math.sqrt(float(np.sum((r * w) ** 2)))
    b = R / (Wh * x ** 0.75); out['okna'].append(dict(x=x, n=int(np.sum(w > 0)), data=Dm, model=M, R=R, b=b, b_iid=eiid / (Wh * x ** 0.75), rel=R / M))
    print('x %.3g  n %6d  Σ L³W %.6g  Σ Q₃W %.6g  R %+.4g (%.2e rel)  b %+.3f ± %.3f (iid)' % (x, np.sum(w > 0), Dm, M, R, R / M, b, eiid / (Wh * x ** 0.75)), flush=True)
bs = np.array([o['b'] for o in out['okna']]); b = float(bs.mean()); sb = float(bs.std(ddof=1) / math.sqrt(len(bs)))
out.update(b=b, sb=sb, z0=b / sb, zres=(b - RES) / sb); out['verdikt'] = 'N' if abs(b / sb) <= 3 else ('O' if abs((b - RES) / sb) <= 3 else 'Z')
print('b %.3f ± %.3f (z od 0 %.1f, od Res %.1f)  signál Res·Ŵ·x^{3/4} u x = 1,25·10⁶: %.1f  verdikt %s' % (b, sb, b / sb, (b - RES) / sb, RES * Wh * 1.25e6 ** 0.75, out['verdikt']))
json.dump(out, open('/home/claude/br79/data/b34.json', 'w'), indent=1)
