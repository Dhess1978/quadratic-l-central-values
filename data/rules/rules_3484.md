# 34.84 (návrh) — nezávislá replikace u 10⁸ a kontrola L(½) u 10¹⁰ (zapsáno 12. 11. 2026 kolem 19:05 UTC, PŘED výpočtem)

David 12. 11. večer: „Mám udělat mezitím nějaký další test nezávisle?“ (Mac je volný, v cloudu běží baterie 34.83).

## 1. Replikace na nových D: nuly D100c (Mac)
- Data: rodina D100c, |D| = 10⁸ + 10⁶ … 10⁸ + 2·10⁶, asi 608 000 L. Nuly počítá ks_nuly_3482.py (SHA-256 2006a966…, divz 32, xmax 4,3).
  L(½) z PARI už jsou ze 34.77 (STRED_D100bc_2026_11_10). Žádné D se nepřekrývá s D100 ani s D100b.
- Postup stejný jako 34.83 (lu8_analyza.py): rozdíl vážené a nevážené hustoty v 5 koších, chyba = větší z (bootstrap 300, promíchání L 200),
  amplituda = průmět na tvar odchylky C100 sudých.
- odv-lu8-rep: ODVOZENO, když u D100c samotné má amplituda sudých a/chyba > 3, (0,80 − a)/chyba > 4 a liché p > 0,001;
  ZAMÍTNUTO, když a/chyba < 2 nebo (0,80 − a)/chyba < 2; jinak NELZE ODVODIT.
- Očekávání zapsané předem: při poklesu |D|^(−θ) s θ = 0,252 (34.83) amplituda asi 0,32 ± 0,07 (stejně jako D100 ∪ D100b, l se liší jen o 0,01);
  při 1/ln asi 0,80.
- Popisně: D100 ∪ D100b ∪ D100c dohromady (asi 1,2 mil. L) a fit bez B10.

## 2. Kontrola hodnot L(½) u 10¹⁰ (Mac, PARI)
- E10 (12 153 L): L(½, χ_D) přes PARI lfun (ks_nuly_3482.py --hodnota-stred) proti vlastní funkcionální rovnici 34.82 (afe_L2.py).
- Bez verdiktu: kontrola v pořádku, když největší relativní rozdíl < 10⁻⁹ u všech D; jinak se E10 přepočítá a 34.82 se znovu vyhodnotí.

Starší verdikty 34.76–34.83 beze změny. Vyhodnocuje se až po doběhnutí celé rodiny.
Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny.

## 3. Ford–Zaharescu: rozdělení {α·γₙ} mod 1 (doplněno 12. 11. 2026 kolem 21:10 UTC, PŘED výpočtem; David: „Souhlasím“)
- Teorie: rovnoměrnost {αγ} mod 1 (Rademacher za RH, Hlawka bezpodmínečně); odchylka druhého řádu jen pro α = ln(p^k)/2π (Ford–Zaharescu 2005,
  Ford–Soundararajan–Zaharescu 2009). Předpověď bez volného parametru z Landauova vzorce (34.43): pro okno nul T₁ < γ ≤ T₂ s N nulami je hustota
  f(u) = 1 − 2w·Σ_{k≥1} Λ(p^k)·p^{−k/2}·cos(2πku), w = (T₂ − T₁)/(2π·N), součet do k = 60.
- Data: prvních 100 000 nul (LMFDB, t ≤ 74 921) jako 1 blok; zeros6 (2 001 052 nul, t ≤ 1,13·10⁶) ve 4 blocích po 500 263; Odlyzkovy bloky 10¹², 10²¹, 10²²
  (po 10 000 nulách, γ = základ + posun, {α·základ} v mpmath na 40 míst) jen popisně (signál prvního harmonického tam ~1 chyba).
- 50 košů na [0, 1). Chyba koše = větší z (√očekávání, směrodatná odchylka 200 nulových modelů). Nulový model: tytéž rozestupy v náhodném pořadí
  (stejná hustota a lokální statistika, bez aritmetiky).
- odv-fz-tvar: α = ln p/2π, p = 2, 3, 5, 7. ODVOZENO, když ve všech 20 kombinacích (5 bloků × 4 prvočísla) p(χ² proti předpovědi) > 0,001 a pro p = 2
  v každém bloku zeros6 p(χ² proti rovnoměrnému) < 10⁻⁶; ZAMÍTNUTO, když aspoň jedna kombinace má p(χ² proti předpovědi) < 10⁻⁴; jinak NELZE ODVODIT.
- odv-fz-kontrola: α = √2/2π, ln 6/2π, ln 10/2π (e^{2πα·k} nikdy mocnina prvočísla) → rovnoměrné rozdělení. ODVOZENO, když ve všech 15 kombinacích
  p(χ² proti rovnoměrnému) > 0,001; ZAMÍTNUTO, když aspoň jedna < 10⁻⁴; jinak NELZE ODVODIT.
- Očekávání zapsané předem: zvlnění prvního harmonického u zeros6 2w·ln p/√p ≈ 8 % (p = 2) až 12 % (p = 7) při šumu koše asi 1 %; u Odlyzkových bloků asi 2 %.
  Popisně: α = ln 4/2π (harmonické 4^k = 2^{2k} s vahou ln 2/2^k), první Fourierův koeficient proti Landauovu m(x) ze 34.43.
