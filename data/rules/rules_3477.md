# 34.77 STŘED V · DOPLŇKY — pravidla verdiktů (zapsáno 10. 11. 2026 večer, před výpočtem)

David 10. 11.: „souhlasím se všemi návrhy“, „zapracuj vše“. Data jako v 34.76 (všechna fundamentální D, 3 ≤ |D| ≤ 10⁷; B10, C100, R10, D100 s prvními nulami).
Polynomy CFKRS s opraveným ocasem (v2). Chyba „z bloků“ = rozptyl součtů po souvislých blocích 20 000; „iid“ = jako u nezávislého výběru.

Poctivě, co už bylo vidět před zápisem: průměr a korelace x₁ s rozvinutím θ_D(γ)/π (studie III) — spočteno 10. 11. večer při kontrole normalizace
(střed x₁ 0,62–0,69; korelace sudé 0,842–0,883, liché 0,833–0,888). H3 a H4 jsou proto dodatečné. Ostatní hodnoty (k = 5, třídy, b₃, zlomkové momenty,
chyby κ₃, malé hodnoty po třídách) zatím spočtené nejsou.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lp-k5 | pátý moment obou znamének sedí s plným polynomem CFKRS | NELZE ODVODIT, když vedoucí koeficient Q₅ z kontury nesedí s a₅·Π j!/(2j)! na 10⁻³; jinak ODVOZENO, když |z_bloky| < 3 celkem (obě znaménka) i v 7 pásmech a |z_iid| < 3 u R10 a D100; jinak ZAMÍTNUTO | jen vedoucí člen |
| odv-lp-x1theta | (dodatečně) korelace ln L s x₁ (rozvinutí θ_D) = USp(2N), N z prvního momentu | ODVOZENO, když |z| < 3 u všech 8 skupin; jinak ZAMÍTNUTO | promíchání x₁ |
| odv-lp-mx1 | (dodatečně) střední první nula (rozvinutí θ_D) = USp(2N), N z prvního momentu | ODVOZENO, když |z| < 3 u všech 8 skupin (chyba dat a DPP); jinak ZAMÍTNUTO | — |
| odv-lp-trida | nadbytek korelace nad USp nesou malá prvočísla: korelace uvnitř tříd (χ_D(2), χ_D(3), χ_D(5)) je blíž USp(2N) než celková | ODVOZENO, když se průměrná |corr − USp| přes 8 skupin zmenší aspoň o polovinu; ZAMÍTNUTO, když se nezmenší; jinak NELZE ODVODIT | celková korelace |
| odv-lp-b3 | sekundární člen b·X^{3/4} u k = 3 (Zhang, Alderson–Rubinstein: b₊ ≈ −0,14, b₋ ≈ −0,07) je při |D| ≤ 10⁷ vidět | odhad b = průměr E₃(X)/X^{3/4} přes X ∈ [2·10⁶, 10⁷] (konce bloků), chyba ze 4 nezávislých oken; ODVOZENO, když |b/σ| > 3 a b je v rámci 3σ od předpovědi; ZAMÍTNUTO, když |b/σ| > 3 a b je mimo; jinak NELZE ODVODIT | b = 0 |
| odv-lp-zlom | hybrid (aritmetika + USp(2N), N z k = 1) popíše momenty s = ½ a s = −½ do 1 % | ODVOZENO, když |data/hybrid − 1| < 0,01 v obou znaménkách v pásmech 10⁶–3·10⁶ a 3·10⁶–10⁷ pro s = ±½; jinak ZAMÍTNUTO | Keating–Snaith s 2N = ln(|D|/π) |
| odv-lp-k3 | třetí kumulant ln L sedí s hybridem v rámci 3 chyb | chyba = větší z iid a jackknife po blocích 20 000; ODVOZENO, když |z| < 3 v obou znaménkách ve třech pásmech |D| ≥ 10⁵; jinak ZAMÍTNUTO | — |
| odv-lp-maletr | schodek malých hodnot proti hybridu je stejný ve všech třídách χ_D(2), χ_D(3) | počty L < 0,03 u |D| ≥ 10⁶ po 9 třídách proti hybridu s podmíněným aritmetickým faktorem, normováno celkovým poměrem; ODVOZENO, když χ² homogenity má p > 0,01; jinak ZAMÍTNUTO | — |

Čeká na data (34.78): šířka bloku u 10⁸ do w = 10⁵ (D100b, D100c z Macu, ks_nuly_3477.py); hustota nul vážená L(½) (nízké nuly B10 … D100 ze studie III).
