# 34.85 STŘED V · VIII (návrh) — pravidla verdiktů (zapsáno 13. 11. 2026 PŘED výpočtem části A i před během S7 na Macu)

David 13. 11.: „využijeme i mac, uděláme to podle tvého návrhu … do max 9 hodin běhu na MAC“.
Proč: jediné nové tvrzení studie V je člen sudých charakterů ve vážené hustotě nízkých nul (amplituda ∝ |D|^−θ, θ = 0,276 ± 0,040).
Článek sám uvádí jeho slabinu: obě implementace Fazzariho předpovědi sdílejí aritmetický faktor, chybu společnou oběma nelze vyloučit,
a amplituda je průmět na tvar viděný v C100. Verze 34.85 to má uzavřet před vydáním studie V: buď člen potvrdit a lokalizovat, nebo ho opravit.
Starší verdikty 34.76–34.84 zůstávají beze změny; výsledek 34.85 se do studie V zapíše takový, jaký vyjde (i kdyby člen padl).

## Část A (cloud, data, která už máme: hodnoty L(½) všech |D| ≤ 10⁷ a u 10⁸, nuly B10, C100, R10, D100, D100b, D100c, E)

Testovací funkce φ(u) = (sin(πu/2)/(πu/2))⁴, u = γ·ln(|D|/π)/(2π); její Fourierova transformace má nosič |ν| ≤ 1, takže prvočíselná strana
explicitního vzorce končí u n ≤ |D|/π. Δ = (průměr Σ_γ φ(u_γ) vážený L(½, χ_D)) − (nevážený průměr), zvlášť sudé a liché.
Chyby: Poissonův bootstrap přes D (300×) společně pro všechny strany.

| id | hypotéza | pravidlo |
|---|---|---|
| odv-lv-ef | explicitní vzorec (věta, ne domněnka) spojuje nuly a prvočísla v našich datech: Δ ze strany nul = Δ ze strany prvočísel (vážené průměry χ_D(n) z naměřených L(½), Λ(n), archimedovský člen) | C100 obě znaménka (a S7, až bude): ODVOZENO, když se obě strany liší o méně než 3 chyby všude; ZAMÍTNUTO, když kdekoli o více než 5 chyb (pak je problém v datech nebo zpracování a člen sudých se nesmí vydat); jinak NELZE ODVODIT. Ocas nad spočtenými nulami (x > 4,3) se doplní z hustot domněnky poměrů a jeho podíl se uvede |
| odv-lv-nez | nezávislá implementace předpovědi dá táž čísla | nový kód: aritmetický faktor domněnky poměrů jako přímá lokální řada (součet přes mocniny p numericky, ne uzavřený tvar z fazzari.py), derivace podle β komplexním krokem, výměnné členy znovu; předpověď vážené − nevážené hustoty v 5 koších pro B10, C100, R10, D100, D100b, D100c, E a obě znaménka. ODVOZENO, když se od dosavadní předpovědi liší nejvýš o 0,05 chyby dat v každém koši; ZAMÍTNUTO, když kdekoli o více než 0,5 chyby (pak se člen sudých přepočítá s novou předpovědí dřív než cokoli dalšího); jinak NELZE ODVODIT |
| odv-lv-hladka | člen sudých je vidět i hladkou testovací funkcí s nosičem ≤ 1 (nejen v koších) | C100: (Δ_nuly − Δ_předpověď)/chyba > 3 u sudých a < 2 u lichých → ODVOZENO; sudé < 2 → ZAMÍTNUTO (pak člen leží ve frekvencích nad nosičem 1, tj. pochází z n > |D|/π — to je lokalizace, ne vyvrácení výsledku v koších); jinak NELZE ODVODIT |

Popisně (lokalizace, bez verdiktu): prvočíselná strana rozložená na n = 2^k, lichá prvočísla p ≤ (|D|/π)^½, lichá prvočísla nad, mocniny p^k s k ≥ 2,
zvlášť sudé a liché, naměřené proti předpovědi téhož rozkladu z nezávislé implementace (zkroucený první moment ⟨L(½)χ_D(n)⟩ podle receptu);
poměr každé části mezi C100 a S7.

## Část B (Mac, nuly nové rodiny S7)

Data: nuly všech základních D s 7,5·10⁶ < |D| ≤ 10⁷ (obě znaménka, asi 1 519 800 L), ks_nuly.py verze 34.85 (první podoba se samotnou S7 měla SHA-256 92ea2e93…; platná verze viz část C;
výpočet beze změny proti 34.82: divz 32, xmax 4,3), hodnoty L(½) z V10 (PARI, 34.72). Postup jako 34.83/34.84: rozdíl vážené a nevážené hustoty v 5 koších x ∈ [0, 1]
proti domněnce poměrů s vahou L(½) (Fazzari 4.8, kusy 34.79) pro rozdělení ln(|D|/π) rodiny; chyba = větší z (bootstrap 300×, promíchání L 200×); amplituda = vážený průmět na tvar odchylky C100 sudých.

Předpověď zapsaná předem (fit bez B10 ze 34.84, M = c·e^{−θ(l − l_C100)}): l_S7 = 14,836, amplituda sudých **0,608 ± 0,045**; model c·l_C100/l dává 0,639
(u 10⁷ se oba modely neliší — bod S7 testuje tvar zákona a přesnost θ, rozlišení mocnina × 1/ln stojí na bodech u 10⁶ a 10⁸). Čekaná chyba měření asi 0,040, tedy a/chyba asi 15.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lv-s7 | odchylka sudých u 10⁷ odpovídá zákonu |D|^−θ ze 34.84 | ODVOZENO, když a/chyba > 5 a \|a − 0,608\| ≤ 2·√(chyba² + 0,045²); ZAMÍTNUTO, když \|a − 0,608\| > 3·√(…) nebo a/chyba < 2; jinak NELZE ODVODIT. Liché: p > 0,001 | promíchání L mezi D |
| odv-lv-tvar | tvar odchylky sudých je u 10⁷ stejný jako v C100 (jen menší) | S7 sudé minus a·tvar C100 (chyby obou), χ² se 4 stupni volnosti: ODVOZENO, když p > 0,01; ZAMÍTNUTO, když p < 0,001; jinak NELZE ODVODIT | promíchání L |

## Část C (Mac, šestý bod u 10⁹ — doplněno 13. 11. po Davidově prodloužení času na Macu na 1,5 dne, stále PŘED výpočtem)

Data: nuly všech základních D s 10⁹ < |D| ≤ 10⁹ + 8·10⁵ (rodina S9L, asi 486 000 L, obsahuje S9), týž skript a nastavení (běží po S7 v jednom příkazu, sada „uzavreni“).
Hodnoty L(½): vlastní vyhlazená funkcionální rovnice v cloudu (jako E u 10¹⁰, váhy počítané po úsecích s relativní šířkou do 10⁻⁵).
Kontrola hodnot (pravidlo jako u E10 ve 34.84): PARI na Macu spočte L(½) pro 2 000 náhodných D z S9L (rodina S9K, semínko 34851); hodnoty z cloudu se použijí,
když největší relativní rozdíl na všech 2 000 D je pod 10⁻⁹; jinak je část C NELZE ODVODIT.
Skript ks_nuly.py verze 34.85 s rodinami S9L a S9K má SHA-256 1b53fe0e093e4c998e2eda9c63594e297ab08666da4cfec8644bc29497a103e8 (nahrazuje verzi s jedinou rodinou S7, výpočet stejný).

Předpověď zapsaná předem (týž fit): l_S9L = 19,579, amplituda sudých **0,164 ± 0,040**; model c·l_C100/l dává 0,484. Čekaná chyba měření asi 0,085:
člen sám bude vidět jen asi na 2 chyby, ale od hodnoty pro 1/ln|D| se liší asi o 3,8 chyby.

| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lv-s9 | u 10⁹ je odchylka sudých tak malá, jak předpovídá zákon \|D\|^−θ, a ne tak velká jako při 1/ln\|D\| | ODVOZENO, když \|a − 0,164\| ≤ 2·√(chyba² + 0,040²) a (0,484 − a)/chyba > 3; ZAMÍTNUTO, když \|a − 0,164\| > 3·√(…) nebo (0,484 − a)/chyba < 1; jinak NELZE ODVODIT. Liché: p > 0,001 | promíchání L mezi D |

Popisně: fit θ se všemi body (C100, R10, S7, D100 ∪ D100b ∪ D100c, S9L, E; i s B10) — po doběhnutí S9L.
Kdyby běh skončil dřív (limit), použijí se všechny hotové bloky S9L (rozhodnuto teď).

Popisně: nový fit θ s body C100, R10, S7, D100 ∪ D100b ∪ D100c, E (a s B10); explicitní vzorec a rozklad části A i na S7.
Vyhodnocuje se až po doběhnutí běhu (žádné průběžné nahlížení). Kdyby byl běh zastaven limitem 9 h, použijí se všechny hotové bloky (rozhodnuto teď).

Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny.

## Doplněk D (13. 11. 2026, po výsledku odv-lv-nez a PŘED daty S7)

**odv-lv-nez = ZAMÍTNUTO** (dosavadní předpověď): nezávislá implementace se od dosavadní liší u sudých charakterů až o 13,3 chyby dat (C100), u lichých nejvýš o 0,49.
Příčina nalezena: ve fazzari.py (34.79) má člen s oběma výměnami (α → −α, β → −β) Γ-podíl gS(½+α+β) místo gS(½+α)·gS(½+β); po opravě jediného faktoru
se dosavadní kód shoduje s novým na 10⁻⁴. Táž chyba je ve fazzari_tridy.py (proto shoda obou na 10⁻¹³). S opravenou předpovědí sedí sudá i lichá data
ve všech rodinách (χ² na 5 koších: sudé B10 1,0, C100 3,0, R10 1,4, D100 9,2, D100 ∪ D100b 2,9, D100c 0,8, všechny u 10⁸ 2,3, E 1,5; liché 0,3–2,6).
Podle pravidla se člen sudých přepočítává s novou předpovědí dřív než cokoli dalšího: člen navíc v datech není.

Proto, zapsáno před daty S7 (a případně S9L):
- Hypotézy odv-lv-s7, odv-lv-tvar, odv-lv-s9 a odv-lv-hladka se vyhodnotí tak, jak jsou zapsané, ale s opravenou předpovědí; čekané výsledky jsou ZAMÍTNUTO (amplituda ≈ 0), protože jejich předpoklad (člen sudých) padl.
- Nová hypotéza **odv-lv-s7k**: opravená předpověď domněnky poměrů s vahou L(½) popisuje rodinu S7 pro obě znaménka. ODVOZENO, když p(χ², 5 košů) > 0,001 u sudých i lichých;
  ZAMÍTNUTO, když p < 10⁻⁴ u kteréhokoli znaménka; jinak NELZE ODVODIT. Totéž pro S9L (odv-lv-s9k), pokud se spočte.
- Verdikty 34.79–34.84, které stojí na dosavadní předpovědi vážené hustoty (odv-lr-vaha, odv-lu-e10, odv-lu-skal, odv-lu8-e8, odv-lu8-skal, odv-lu8-rep a popisné fity θ),
  se v laboratoři 34.85 přehodnotí s opravenou předpovědí; původní verdikty zůstanou zapsané s poznámkou „chyba implementace, opraveno 34.85“.

## Doplněk E (14. 11. 2026, PŘED výpočtem nul S9L a před hodnotami L(½) u 10⁹)

David 14. 11.: „když ho spočítáme, budeme tam mít data do budoucna“ → rodina S9L se počítá (pro úplnost řady vodičů 10⁴ … 10¹⁰ a jako data pro další práci);
vydání studie V počká. Po opravě předpovědi (doplněk D) už bod u 10⁹ netestuje žádný zákon škálování.

- Data: beze změny podle části C — nuly všech základních D s 10⁹ < |D| ≤ 10⁹ + 8·10⁵ (S9L), skript ks_nuly_3485.py (SHA-256 1b53fe0e…), xmax 4,3, divz 32.
  Kdyby byl běh přerušen, použijí se všechny hotové bloky (jak je rozhodnuto v části C); vyhodnocuje se až po skončení běhu, bez průběžného nahlížení.
- Hodnoty L(½): vlastní vyhlazená funkcionální rovnice v cloudu, váhy po úsecích |D| s relativní šířkou do 10⁻⁵; použijí se, když se na všech 2 000 D
  kontrolního vzorku S9K (PARI na Macu) liší relativně o méně než 10⁻⁹; jinak je bod u 10⁹ NELZE ODVODIT (pravidlo části C).
- Vyhodnocení: týž postup a týž skript jako u S7 (s7_analyza.py: 5 košů na x ∈ [0, 1], chyba = větší z bootstrapu 300× a promíchání L 200×, amplituda na tvaru
  dřívější odchylky C100, předpověď z nezávislé implementace se 40 koši l).
- Hypotéza **odv-lv-s9k** (z doplňku D): opravená předpověď popíše S9L pro obě znaménka — ODVOZENO, když p(χ², 5 košů) > 0,001 u sudých i lichých;
  ZAMÍTNUTO, když p < 10⁻⁴ u kteréhokoli znaménka; jinak NELZE ODVODIT.
- Hypotéza odv-lv-s9 (část C) se vyhodnotí tak, jak je zapsaná, s opravenou předpovědí; čekaný výsledek ZAMÍTNUTO (předpoklad — člen sudých — padl).
  Pro záznam: chybná předpověď by podle zákona ze 34.84 dala amplitudu 0,164 ± 0,040; čekaná chyba měření asi 0,085, s opravenou předpovědí čekáme 0 v chybě.
- Explicitní vzorec se na S9L NEVYHODNOCUJE (rozhodnuto teď): strana prvočísel by šla do n ≤ |D|/π ≈ 3,2·10⁸, o několik řádů dráž než u S7.
  Hypotéza odv-lv-ef zůstává vyhodnocena na C100 a S7, jak ji pravidlo části A jmenuje.
- Popisně (bez verdiktu): střední první nula a ⟨L(½)⟩ proti předpovědi, χ² se samotným bootstrapem, počet různých D s nulami po přidání S9L.
- Výsledek vstoupí do laboratoře jako verze 34.85b a do studie V (další řádek tabulky 1, bod v obrázcích 6 a 7, věta v oddílu 4.6).
- Upřesnění k odv-lv-s9 (14. 11., stále PŘED nulami S9L; hodnoty L(½) u 10⁹ se v cloudu právě počítají a nebyly nijak vyhodnoceny): věta výše „čekaný výsledek ZAMÍTNUTO“ neplatí obecně. Při čekané chybě
  amplitudy asi 0,085 pravidlo části C nerozliší amplitudu 0 od 0,164 (rozdíl 1,9 chyby): vyjde-li a ≈ 0, je |a − 0,164| ≤ 2·√(chyba² + 0,040²) ≈ 0,19 i
  (0,484 − a)/chyba > 3 splněno a pravidlo dá formálně ODVOZENO, přestože jeho předpoklad (člen sudých) padl. Výsledek odv-lv-s9 se proto uvede jen pro záznam
  (číslo a formální výrok pravidla) a NEVYKLÁDÁ se; o bodu u 10⁹ rozhoduje jen odv-lv-s9k. Rozlišení mezi opravenou a chybnou předpovědí se uvede popisně
  jako (amplituda s chybnou předpovědí − amplituda s opravenou)/chyba.
- Stav 14. 11. 08:53 UTC (před nulami S9L): kontrola hodnot L(½) proti vzorku S9K prošla — největší relativní rozdíl 8,4·10⁻¹³ na všech 2 000 D (mez 10⁻⁹), hodnoty
  všech 486 344 základních D rodiny S9L jsou tím uvolněné. Skripty vyhodnocení jsou zmrazené před příchodem nul: s9l_analyza.py SHA-256 5d42f0bb014f32c1…c4da01f,
  s9l_sha.py 9509d786836c02a8…861c5046, s9k_kontrola.py f4eef7ae487a5a95…45d884, afe_L3.py b7174b8573f1d847…874776, s9l_D.py db5cbeac7ece73db…0ff5d7c.
- Popisně navíc (zapsáno teď, před výpočtem; bez pravidla a bez verdiktu; na nulách nezávisí): momenty k = 1 … 5 hodnot L(½) rodiny S9L proti úplným polynomům
  CFKRS týmž postupem jako H0 u 10⁸ ve 34.78 (zbytek L^k − Q_k(ln|D|), chyba ze součtů po blocích šířky 2 000 v |D|, pro srovnání jen vedoucí člen) — skript s9l_momenty.py.
  Do studie V vstoupí jako popisná věta, ne jako tvrzení pod pravidlem.
