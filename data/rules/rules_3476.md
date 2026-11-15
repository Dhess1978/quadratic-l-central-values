# 34.76 STŘED V — pravidla verdiktů (zapsáno 10. 11. 2026 před sestavením listu)

Poctivě: pravidla jsou zapsaná **po předběžném zpracování v cloudu** (10. 11., STAV_STUDIE_V_ZPRACOVANI.md), ne naslepo před výpočtem.
Většina měření byla v zadání (ZADANI_MOMENTY_STUDIE_V.md, 8. 11.) jen slovně; prahy doplňujeme teď a list je přepočítá živě.
Hypotézy H5 (šířka bloku) a H7 (exponent malých hodnot) jsou nové z 10. 11. odpoledne — zapsány po prvním pohledu na data, proto
jsou v textu listu označené jako dodatečné (post hoc).

Data: V5 ∪ V10 = všechny fundamentální diskriminanty 3 ≤ |D| ≤ 10⁷ (6 079 285 L), rodiny B10, C100, R10, D100 se studie III.
Chyba „empirická“ = z rozptylu součtů reziduí po souvislých blocích 20 000 v |D|; „iid“ = jako u nezávislého výběru.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lh-cfkrs | momenty k = 1–4 obou znamének sedí s plným polynomem CFKRS | ODVOZENO, když |z_emp| < 3 u všech 8 celkových součtů, |z_emp| < 3 u všech 56 pásem (7 pásem |D| × 2 znaménka × 4 k), u R10 a D100 |z_iid| < 3 (16 hodnot); jinak ZAMÍTNUTO | jen vedoucí člen (H2) |
| odv-lh-ks | samotný vedoucí člen Keatinga–Snaithové popisuje momenty do 10 % | ODVOZENO, když |data/KS − 1| < 0,1 u všech 8 celkových průměrů; jinak ZAMÍTNUTO | plný polynom CFKRS |
| odv-lh-parita | rozdíl průměrů L lichých a sudých = A₀·π/2 (Γ-faktory) | ODVOZENO, když |z| < 3 s empirickou chybou; jinak ZAMÍTNUTO | rozdíl 0 |
| odv-lh-bloky | souvislé bloky 20 000 kolísají jako náhodné výběry téže velikosti | ODVOZENO, když poměr sd bloků / sd náhodných výběrů leží v [0,8; 1,25] u všech 8 (znaménko × k); jinak ZAMÍTNUTO | náhodné podmnožiny (400×), R10 |
| odv-lh-sirka | útlum bloků působí periodicita χ_D(n) v D: s šířkou bloku w roste (poměr klesá) a u větších |D| je při téže w slabší | ODVOZENO, když poměr klesá monotonně s w (šířky s ≥ 20 bloky) v obou znaménkách pro k = 1, 2 u V (5·10⁶ < |D| ≤ 10⁷) i D100, a u D100 je poměr větší než u V při w = 500, 1 000, 2 000 ve všech 12 srovnáních; NELZE ODVODIT, když platí jen první část; jinak ZAMÍTNUTO | iid poměr 1 při každé w |
| odv-lh-rozptyl | rozptyl ln L sedí s hybridem (aritmetika + USp(2N), N z prvního momentu, bez volného parametru) | ODVOZENO, když |var_data/var_hybrid − 1| < 0,01 v obou znaménkách ve třech pásmech |D| ≥ 10⁵; jinak ZAMÍTNUTO | Keating–Snaith s 2N = ln(|D|/π) |
| odv-lh-male | P(L < ε) ∝ ε^{3/2} pro malá ε (pól E[Λ^s] v s = −3/2 u USp) | ODVOZENO, když odhad exponentu (max. věrohodnost, ε₀ = 0,03) je u sudých, lichých i všech v rámci 3 chyb od 3/2; jinak ZAMÍTNUTO | bez odpuzování: exponent 1 |
| odv-lh-x1 | korelace ln L s první nulou x₁ = USp(2N) s N z prvního momentu | ODVOZENO, když |corr_data − corr_USp(N_k1)| < 3 chyby u všech 8 skupin (4 rodiny × znaménko; chyba = data a simulace DPP dohromady); jinak ZAMÍTNUTO | promíchání x₁ (corr 0) |
| odv-lh-nevym | žádné L(½, χ_D) = 0 pro 3 ≤ |D| ≤ 10⁷ | ODVOZENO, když nejmenší hodnota > 0 a počet ≤ 0 je 0; jinak ZAMÍTNUTO | — |

Kontroly bez verdiktu: SHA-256 bloků, síto fundamentálních diskriminantů živě, L(½) živě (vyhlazená AFE) proti PARI, mpmath 30 cifer,
Q₁ analyticky živě proti konturovému integrálu, vedoucí koeficienty a_k·Π j!/(2j)! živě, USp(2N) přesné momenty živě, DPP živě.

Vysvětlení dřívější poznámky: „hybrid s nezávislým aritmetickým faktorem dává 0,76–0,81 → zamítnuto“ (STAV 10. 11. dopoledne) se
**nepřebírá** — výpočet se nepodařilo zopakovat a aritmetický „rozptyl“ κ₂ = −0,289 je záporný, takže nezávislá náhodná veličina
s tímto faktorem neexistuje. Místo toho H8 (odv-lh-x1) s jasným pravidlem.
