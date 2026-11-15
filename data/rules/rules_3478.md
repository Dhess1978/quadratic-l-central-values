# 34.78 STŘED V · III — pravidla verdiktů (zapsáno 11. 11. 2026 odpoledne, PŘED výpočtem H1–H4)

David 11. 11.: „Pokračuj co potřebuješ“ (po návrhu tří doplňků ke studii V: vážená hustota nul, rozklad P_X · Z_X, škálování útlumu bloků) a dodání PDF Aldersona a Rubinsteina.
Data: L(½, χ_D) všech fundamentálních D, 3 ≤ |D| ≤ 10⁷ (V5, V10); B10, C100, R10, D100 s nízkými nulami (studie III, cloud br63/KS = Mac, manifest MD5 4e0a545c… a ccc49e9f…; nuly do x = 4,3, krok lfunzeros střední rozestup/8 — známá neúplnost blízkých dvojic, odv-kh-uplnost); nově D100b + D100c (|D| = 10⁸ + 10⁵ … 10⁸ + 2·10⁶, Mac, běh dlouhy8).
Poctivě, co je už spočtené: srovnání našich Q₁–Q₅ s tabulkami Aldersona a Rubinsteina (11. 11. před zápisem, claude/AR_SROVNANI_3478.md) — není hypotéza, jde do listu jako SHODA S ČÍSLEM ZVENKU. Nic z H1–H4 spočtené není.
Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny.

## Ověření dat z Macu (kontrola, ne hypotéza; při selhání se nic dál nepočítá a David dostane zprávu)
SHA-256 všech 950 bloků proti manifestu; vlastní síto fundamentálních D v rozsahu = množina D v datech (žádné chybějící, navíc ani dvakrát; navázání na D100 u 10⁸ + 10⁵); 2 000 náhodných + 20 nejmenších hodnot vlastní AFE (bez PARI) do 10⁻¹² relativně; 200 hodnot mpmath na 30 cifer.

## Hypotézy
| id | hypotéza | pravidlo | nulový model |
|---|---|---|---|
| odv-lq-cfkrs8 | momenty k = 1–5 u D100b + D100c (asi 1,15 milionu L, |D| ≈ 10⁸) sedí s plným polynomem CFKRS | z = (data − CFKRS)/chyba z bloků 2 000; ODVOZENO, když |z| < 3 u všech k a obou znamének; ZAMÍTNUTO, když někde |z| ≥ 4; jinak NELZE ODVODIT | jen vedoucí člen |
| odv-lq-skal | útlum souvislých bloků závisí na |D| jen přes u = w/√|D| | ρ(w) = sd normovaných součtů bloků (L^k − Q_k, k = 1, 2) / iid; pásmo A |D| ∈ [10⁶, 3·10⁶] (V), pásmo B D100b + D100c; w_B ∈ {10³, 3·10³, 10⁴, 3·10⁴}, w_A = w_B·√(2·10⁶/1,05·10⁸); chyba sd z počtu bloků (1/√(2(n−1)) relativně); 16 srovnání (4 šířky × 2 znaménka × k = 1, 2); ODVOZENO, když max |z| < 3; ZAMÍTNUTO, když max |z| ≥ 4; jinak NELZE ODVODIT. Popisně: w = 10⁵ (19 bloků), srovnání se škálováním „w samo“ a „w/|D|“ a s logaritmickou variantou z mechanismu | náhodné výběry téže velikosti (ρ = 1) |
| odv-lq-vaha | posun střední první nuly při vážení L(½)^k (k = 1, 2) odpovídá USp(2N) váženému Λ(1)^k | Δ_k = E[L^k x₁]/E[L^k] − E[x₁], x₁ rozvinutím θ_D (studie III); N_x1 ze shody neváženého E[x₁] (interpolace mezi celými N, DPP 200 000 matic na N); 8 skupin (B10, C100, R10, D100 × znaménko); chyba dat bootstrapem přes D (1 000×), chyba RMT ze vzorku; ODVOZENO, když max |z| < 3 (16 hodnot); ZAMÍTNUTO, když max |z| ≥ 4; jinak NELZE ODVODIT | promíchání L mezi D uvnitř skupiny (Δ = 0) |
| odv-lq-hust | rozdíl hustoty nízkých nul vážené L(½) a nevážené odpovídá USp(2N_x1) (Fazzari: v limitě W¹_USp − W_USp) | rozvinuté nuly x ∈ [0, 3], 15 košů po 0,2; rozdíl vážené − nevážené hustoty, data proti USp(2N_x1) vážené Λ(1); χ² s chybami z bootstrapu přes D; ODVOZENO, když p > 0,001 ve všech 8 skupinách; ZAMÍTNUTO, když p < 10⁻⁴ aspoň ve 2 skupinách; jinak NELZE ODVODIT. Popisně: limita W¹_USp − W_USp | promíchání L mezi D |
| odv-lq-split | hybridní Eulerův–Hadamardův součin: momenty se dělí, E[L^k] = E[P_X^k]·E[Z_X^k] | P_X = Π_{p ≤ X, p ∤ D}(1 − χ_D(p)/√p)⁻¹, Z_X = L/P_X; |D| ∈ [3·10⁶, 10⁷], obě znaménka; X ∈ {10, 30, 100, 300}; poměr S_k(X) = E[L^k]/(E[P_X^k]E[Z_X^k]) pro k = 1, 2, 3, chyba z bloků 20 000; ODVOZENO, když |S_k − 1| < 0,01 pro všechna k a X ≤ 100 v obou znaméncích; ZAMÍTNUTO, když |S_k − 1| > 0,03 pro některé X ≤ 100; jinak NELZE ODVODIT. Popisně: korelace ln P_X a ln Z_X, X = 300. Očekávání zapsané předem: u ζ laboratoř dělení na dosažitelných výškách zamítla (odv-hs-split, poměr 0,80 už u X = 10) | — |

## Čísla zvenku
- Alderson–Rubinstein 2012, tabulky 1–2: naše Q₁–Q₅ = jejich předpovědi do 6,5·10⁻⁸ (X do 5·10¹⁰); b₊, b₋ (Zhang) jen popisně.
- Fazzari (arXiv 2109.07244): W¹_USp(x) = 1 + sin(2πx)/(2πx) − 2 sin²(πx)/(πx)² (podmíněně RH a domněnkou poměrů).
- Gonek–Hughes–Keating, Bui–Florea: rozklad L = P_X·Z_X a domněnka dělení (asymptoticky).
