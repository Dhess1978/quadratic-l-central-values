# 34.83 STŘED V · VI (návrh) — pravidla verdiktů (zapsáno 12. 11. 2026 kolem 10:00 UTC, PŘED výpočtem nul D100b)

David 12. 11. (při běhu ověření 34.82): nabídl volný výkon Macu pro kontrolní běh („když mi dáš PY soubor, mohu ti udělat kontrolní běh“).
Proč: ve 34.82 stojí zamítnutí poklesu 1/ln|D| (odv-lu-skal) hlavně na rodinách B10 a C100. Bez B10: θ = 0,22 ± 0,10 a model 1/ln|D| má p 0,05.
Nejsilnější nezávislý bod je u 10⁸: rodina D100 má nuly jen pro 60 800 L, ale hodnoty L(½) z PARI už jsou i pro D100b (34.77, 547 000 L).
Chybí jen nuly D100b.

Data (až doběhnou):
- Nuly D100b: Davidův Mac, ks_nuly.py verze 34.82 (SHA-256 2006a966a50180a67875375c32c7876821c84c1908746eaa556f8f82ee2f9554, beze změny výpočtu
  proti 34.77), divz 32, xmax 4,3. Úplnost se ověří přes ks_kontrola.py a SHA bloků.
- Nuly D100 ze sady velky; L(½) z PARI lfun (34.72 pro D100, 34.77 pro D100b).
- Spojená skupina D100 ∪ D100b: |D| = 10⁸ … 10⁸ + 10⁶, asi 608 000 L, obě znaménka zvlášť.

Postup stejný jako ve 34.82: rozdíl vážené a nevážené hustoty v 5 koších x ∈ [0, 1] proti domněnce poměrů s vahou L(½) (Fazzari, rov. 4.8, kusy 34.79)
pro rozdělení ln(|D|/π) skupiny. Chyba = větší z (Poissonův bootstrap 300×, rozptyl promíchání L 200×). Amplituda = vážený průmět na tvar odchylky
C100 sudých.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lu8-skal | odchylka sudých klesá s vodičem jako mocnina \|D\|, ne jako 1/ln\|D\| — i bez B10 | body C100, R10, D100 ∪ D100b, E (sudé; B10 vynechán); model M = c·e^{−θ(l − l_C100)}, model L = c·l_C100/l; ODVOZENO, když M má p > 0,01 a L p < 0,001; ZAMÍTNUTO, když L má p > 0,01; jinak NELZE ODVODIT. Očekávání zapsané předem: platí-li mocnina s θ = 0,245, vyjde amplituda u 10⁸ asi 0,35 ± 0,07 a L p kolem 2·10⁻⁵; platí-li 1/ln, amplituda asi 0,80 a L p kolem 0,08 | a = 0 |
| odv-lu8-e8 | u 10⁸ odchylka sudých ještě trvá (je to nižší člen, ne šum), liché sedí | D100 ∪ D100b; ODVOZENO, když amplituda sudých a/chyba > 3 a liché mají p > 0,001; ZAMÍTNUTO, když a/chyba < 2; jinak NELZE ODVODIT. Očekávání: při mocnině a/chyba asi 4,8 | promíchání L mezi D |

Popisně: fit se všemi pěti body (včetně B10), liché amplitudy, θ ± chyba.
Vyhodnocuje se až celá D100b (žádné průběžné nahlížení a zastavování podle výsledku). Starší verdikty 34.76–34.82 beze změny.
Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny.
