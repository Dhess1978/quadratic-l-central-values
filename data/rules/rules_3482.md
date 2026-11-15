# 34.82 STŘED V · V — pravidla verdiktů (zapsáno 12. 11. 2026 v 08:00 UTC, PŘED výpočtem)

David 12. 11.: „pokračuj ve studii V“ a volba „Nový list + sloučení“ (34.82 = 34.81 + STŘED V · IV + nový list k otevřené otázce).
Otevřená otázka z 34.79: vážená hustota nízkých nul (Fazzari, domněnka poměrů s vahou L(½)) sedí u lichých D, u sudých D se odchyluje
(B10 a C100 jistě, R10 a D100 v rozptylu nulového modelu), ve všech třídách χ_D(2). Literatura: Fazzari (arXiv 2109.07244), Berry–Keating 2002
(J. Phys. A 35 L1, záporné momenty a shluky téměř degenerovaných hladin), Keating–Snaith 2000.
Nová data: L(½, χ_D) pro všechna D rodin E10 a E100 (|D| od 10¹⁰ do 10¹⁰ + 2·10⁵, nízké nuly z balíků KS_NULY), počítané v cloudu vlastní
vyhlazenou funkcionální rovnicí (2 Σ χ_D(n) n^{−½} Γ(c, πn²/|D|)/Γ(c), c = ¼ sudé, ¾ liché, n ≤ √(55|D|/π) + 2), kontrola PARI lfun na náhodném vzorku.
Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny; staré verdikty beze změny.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lu-e10 | vážená hustota nízkých nul u |D| ≈ 10¹⁰ odpovídá domněnce poměrů s vahou L(½) (Fazzari, rov. 4.8, nižší členy, kusy z 34.79) u sudých i lichých D | E10 + E100, sudé a liché zvlášť; rozdíl vážené − nevážené hustoty v 5 koších po 0,2 v x ∈ [0, 1] proti předpovědi pro rozdělení ln(|D|/π) skupiny; chyba = větší z (Poissonův bootstrap 300×, rozptyl promíchání L 200×); ODVOZENO, když p > 0,001 u obou znamének; ZAMÍTNUTO, když p < 10⁻⁴ aspoň u jednoho; jinak NELZE ODVODIT. Očekávání zapsané předem: kdyby odchylka sudých klesala jako 1/ln|D| od C100, čeká se u sudých χ² kolem 30 (ZAMÍTNUTO); kdyby jako mocnina |D|, χ² kolem 5 (ODVOZENO) | promíchání L mezi D (rozdíl ≈ 0) |
| odv-lu-skal | odchylka sudých od předpovědi klesá s vodičem jako mocnina |D|, ne jako 1/ln|D| | amplituda a_g = vážený průmět (data − předpověď) na tvar odchylky C100 sudých (váhy 1/chyba²), chyba z chyb košů; 5 bodů B10, C100, R10, D100, E (sudé) proti l = střední ln(|D|/π); model M = c·e^{−θ(l − l_C100)} (2 parametry) a model L = c·l_C100/l (1 parametr); ODVOZENO, když M má p > 0,01 a L p < 0,001; ZAMÍTNUTO, když L má p > 0,01; jinak NELZE ODVODIT. Popisně: θ ± chyba, liché amplitudy týmž postupem | a = 0 |
| odv-lu-male | konstantu malých hodnot P(L < ε) ≈ C·ε^{3/2} určuje shluk u středu (nejnižší nula a její zrcadlo; Berry–Keating 2002): hybrid a(−3/2)·(2/3)·Res_{s=−3/2} E[Λ^s] s N z rozdělení první nuly (N_x1) | skupiny B10, C100, R10, D100, E × znaménko (10 skupin); ε = 0,03, počet hodnot pod ε proti předpovědi C·ε^{3/2}·n; N_x1 ze střední x₁ (USp(2N), Jacobiho tabulka 34.78, lineárně mezi celými N; u E stejně z jejích nul); a(−3/2) = Eulerův součin rodiny (34.76); chyba Poissonova z předpovědi; ODVOZENO, když |z| < 3 ve všech 10; ZAMÍTNUTO, když aspoň jedno |z| ≥ 4; jinak NELZE ODVODIT. Srovnání: N z prvního momentu (34.76: hybrid o 12–20 % víc malých hodnot). Popisně ε = 0,01 a 0,1 | N z prvního momentu (34.76) |

Kontroly bez verdiktu: L z vlastní funkcionální rovnice proti PARI lfun (náhodný vzorek, rozdíl do 10⁻⁹ relativně), úplnost E10 + E100 (všechna D z balíků nul mají L), žádné L ≤ 0;
Fazzariho předpověď u l ≈ 22 z kusů 34.79 (kontrola: nevážená část = vzorec 34.62).
