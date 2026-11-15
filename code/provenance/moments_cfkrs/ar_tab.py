# Tabulky 1 a 2 z Alderson–Rubinstein (arXiv 1110.0253v2), přepsáno z pdftotext -layout; sloupce: k, součet dat, předpověď CFKRS (jejich), R, Δ
import re, json
T = open('/home/claude/ar_pdf/ar.txt', encoding='utf-8').read().splitlines()
def tab(a, b):
    out = []
    for l in T[a:b]:
        p = l.split()
        if len(p) == 5 and p[0].isdigit() and 1 <= int(p[0]) <= 8:
            f = lambda s: float(s[:-1] if s.endswith('.') else s)
            out.append(dict(k=int(p[0]), suma=f(p[1]), pred=f(p[2]), R=f(p[3]), D=f(p[4]), s_suma=p[1], s_pred=p[2]))
    return out
i1 = next(i for i, l in enumerate(T) if l.strip().startswith('2.1') and 'Tables and Figures' in l); i2 = next(i for i, l in enumerate(T) if l.startswith('Table 1:')); i3 = next(i for i, l in enumerate(T) if l.startswith('Table 2:'))
A = dict(minus=tab(i1, i2), plus=tab(i2, i3))
for s, X0 in (('minus', [1e10 * j for j in range(1, 6)]), ('plus', [2e9 * j for j in range(1, 6)])):
    for r in A[s]: pass
    for k in range(1, 9):
        rows = [r for r in A[s] if r['k'] == k]; assert len(rows) == 5, (s, k, len(rows))
        for r, X in zip(rows, X0): r['X'] = X
print({s: len(v) for s, v in A.items()})
# kontrola přepisu: R = suma/pred a Δ = suma − pred
bad = [(s, r['k'], r['X']) for s in A for r in A[s] if abs(r['suma'] / r['pred'] - r['R']) > 2e-9 * max(1, abs(r['R'])) or abs((r['suma'] - r['pred']) - r['D']) > 1e-6 * abs(r['suma']) + 1]
print('nesedí vnitřní kontrola:', bad)
json.dump(A, open('ar_tabulky.json', 'w'), indent=1)
