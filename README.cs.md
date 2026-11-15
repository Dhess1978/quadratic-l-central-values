# Hodnoty kvadratických Dirichletových L-funkcí ve středu pro všechny fundamentální diskriminanty do 10⁷ a u 10⁸, 10⁹ a 10¹⁰ (studie V)

**David Hess (MAYAN ALFA V1008 LAB v.34.85b) · 2026** · [English version](README.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24256378.svg)](https://doi.org/10.5281/zenodo.24256378) Zenodo: [doi:10.5281/zenodo.24256378](https://doi.org/10.5281/zenodo.24256378)

Pátá studie řady, po studiích *Berry saturation and the statistics of the Riemann zeros* (studie I, [doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328)), *Beyond the diagonal approximation* (studie II, [doi:10.5281/zenodo.23928017](https://doi.org/10.5281/zenodo.23928017)), *Low-lying zeros of 795,381 quadratic Dirichlet L-functions* (studie III, [doi:10.5281/zenodo.23991165](https://doi.org/10.5281/zenodo.23991165)) a *Three-point correlation of the Riemann zeros at finite height* (studie IV, [doi:10.5281/zenodo.24094318](https://doi.org/10.5281/zenodo.24094318)). Měří hodnoty L(½, χ_D) kvadratických Dirichletových L-funkcí ve středu a spojuje je s nízkými nulami týchž L-funkcí (studie III). Data: všechny fundamentální diskriminanty 3 ≤ |D| ≤ 10⁷ (6 079 285 hodnot, PARI/GP), všech 1 215 839 diskriminantů s 10⁸ < |D| ≤ 10⁸ + 2·10⁶ (PARI/GP), všech 486 344 s 10⁹ < |D| ≤ 10⁹ + 8·10⁵ a všech 121 557 s 10¹⁰ < |D| ≤ 10¹⁰ + 2·10⁵ (vlastní vyhlazená funkcionální rovnice, ověřená v PARI na 2 000 náhodných diskriminantech u 10⁹ a na všech 12 153 diskriminantech prvního pásma u 10¹⁰); pro 4 042 611 z nich i nízké nuly (nuly 3 161 211 L-funkcí byly spočteny pro tuto studii). Většina tvrzení byla testována podle pravidel zapsaných před výpočtem (`data/rules/`).

| výsledek | hodnota |
|---|---|
| momenty k = 1–5 proti plným polynomům CFKRS, 7 pásem \|D\| ≤ 10⁷ × 2 znaménka × 5 (70 kombinací) | největší odchylka 0,65 směrodatné chyby; celkem 6·10⁻⁷ – 3·10⁻⁴ relativně |
| totéž u 10⁸ (D100b ∪ D100c, 1 155 037 hodnot) | největší odchylka 0,57 chyby |
| popisně (bez zapsaného pravidla) totéž u 10⁹ (S9L, 486 344 hodnot) | největší odchylka 1,15 chyby |
| samotný vedoucí člen (Keating–Snaithová s aritmetickým faktorem) | mýlí se o 9 % (k = 1) až o faktor 6·10⁴ (k = 5) |
| naše polynomy CFKRS proti předpovězeným hodnotám Aldersona a Rubinsteina (k = 1–5, obě znaménka, do X = 5·10¹⁰) | 6,5·10⁻⁸ |
| rozdíl prvních momentů lichých − sudých: data / A₀π/2 | 1,10649 ± 0,00060 / 1,10654 |
| součty přes bloky 20 000 po sobě jdoucích \|D\| proti náhodným výběrům stejné velikosti (k = 1–4) | rozptyl 0,11–0,33 rozptylu nezávislých hodnot (náhodné výběry 0,94–1,04); klesá se šířkou w; škálování s w/√\|D\| i s ln w/ln √\|D\| zamítnuto |
| rozptyl ln L proti hybridu (aritmetický faktor × USp(2N), N z prvního momentu, bez volného parametru) | do 0,26 % |
| exponent malých hodnot, P(L < ε) ∝ ε^β | β = 1,518 ± 0,015 (USp: 3/2); o 11–17 % méně malých hodnot než hybrid |
| nejmenší hodnota; hodnoty ≤ 0 | 6,04·10⁻⁵ u D = −115 147; žádné |
| korelace ln L(½) s první nulou: data / USp(2N) s N z prvního momentu | 0,83–0,89 / 0,81–0,89 (data výš v 7 z 8 skupin, až o 17 chyb) |
| N_ef z domněnky poměrů: střední první nula / její posun při vážení L(½) | vysvětlena / stále o 9–16 % nad USp(2N_ef) |
| vážená jednohladinová hustota (domněnka poměrů s vahou L(½), Fazzari), obě znaménka, 10⁴–10¹⁰ | χ² 0,2–4,3 na 5 koších ve 20 skupinách (nejmenší p 0,51; chyby jsou opatrné, oddíl 4.6) |
| rodina spočtená až po opravě podle pravidla zapsaného před jejími daty (S7: všechna 7,5·10⁶ < \|D\| ≤ 10⁷, 1 519 830 L, 6 592 004 nul) | χ² 0,21 / 0,64 (p 0,999 / 0,986); amplituda sudých 0,015 ± 0,040 (chybná implementace by dala 0,62, jak předpověděl její zákon škálování: 0,608 ± 0,045) |
| rodina u 10⁹ spočtená až po opravě podle pravidla zapsaného před jejími daty (S9L: všechna 10⁹ < \|D\| ≤ 10⁹ + 8·10⁵, 486 344 L, 2 099 543 nul; hodnoty vlastní funkcionální rovnicí) | χ² 1,75 / 1,31 (p 0,883 / 0,934); amplituda sudých −0,088 ± 0,084 (chybná implementace by dala 0,13) |
| oprava před vydáním: naše první implementace této předpovědi (jeden chybný Γ-faktor) vytvářela zdánlivý člen sudých charakterů klesající jako \|D\|^(−1/4); nezávislá implementace podle pravidla zapsaného předem chybu našla | rozdíl až 13,3 chyby před opravou, 0,005 chyby po ní; amplitudy sudých nyní 0,065 ± 0,059 (C100), −0,038 ± 0,050 (10⁸) |
| explicitní vzorec na C100 (607 941 L) a S7 (1 519 830 L): strana nul proti straně prvočísel z naměřených L(½), vážené − nevážené | shoda na 2,2 / 1,6 chyby na C100 a 2,0 / 0,1 na S7 (sudé / liché) |
| L(½) u 10¹⁰: vlastní funkcionální rovnice proti PARI, 12 153 hodnot | největší relativní rozdíl 4,3·10⁻¹¹ |
| příloha: Ford–Zaharescu na 2 001 052 nulách ζ, {αγ} mod 1 proti Landauovu vzorci | nejmenší p 0,66 (α = ln p/2π), rovnoměrnost pro p = 2 vyloučena; kontroly rovnoměrné |
| nezávislá kontrola (`code/independent_check_v.py`, bez sdíleného kódu) | 12 z 12 kontrol prošlo (`code/independent_check_v.log`) |

Polynomy CFKRS, exponent 3/2, hybridní modely a vážená hustota domněnky poměrů jsou známá teorie (domněnky, kromě dokázaných částí); příspěvkem jsou společná data hodnot a nul, jev souvislých rodin u momentů, srovnání bez volných parametrů a test vážené hustoty domněnky poměrů pro obě znaménka na 4,0 milionu L-funkcí, včetně dvou rodin (u 10⁷ a u 10⁹) spočtených až po opravě podle pravidel zapsaných před jejich daty, s kontrolou explicitním vzorcem na dvou rodinách. Zdánlivý člen nižšího řádu u sudých charakterů, nalezený s naší první implementací předpovědi, byl chybou této implementace; popisuje ho oddíl 4.6 článku. Shoda je numerická podpora domněnky, ne důkaz; nic zde se netýká zobecněné Riemannovy hypotézy.

## Obsah

| složka | co |
|---|---|
| `paper/` | článek `central_values_quadratic_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); česká verze v `paper/cs/` |
| `figures/` | obrázky 1–9 anglicky (`_en`) a česky (`_cz`), PNG 300 dpi a PDF |
| `data/central_values/` | L(½, χ_D) pro všechny diskriminanty tabulky 1 článku (NumPy `.npz`, pole `D` a `L` seřazená podle \|D\|; části pod 20 MB), `index.json` s počty a SHA-256 |
| `data/results/` | `study_v_results.json`: všechna čísla článku (momenty po pásmech, bloky, ln L, malé hodnoty, první nula, hybridní součiny, vážená hustota s opravenou předpovědí, rodiny S7 a S9L, explicitní vzorec na C100 a S7, dřívější čísla chybné implementace pro záznam, Ford–Zaharescu, verdikty) |
| `data/cfkrs/` | koeficienty polynomů CFKRS Q_k(ln\|D\|), k = 1–5, obě znaménka; tabulky Aldersona a Rubinsteina použité ke srovnání |
| `data/manifests/` | manifesty běhů (SHA-256 každého bloku hodnot a nul, prostředí, otisk skriptu) |
| `data/rules/` | pravidla zapsaná před výpočtem, s daty (česky, jak byla zapsána) |
| `code/` | `figures_v.py` (obrázky ze souboru výsledků), `independent_check_v.py` a jeho log, `weighted_density_independent.py` (nezávislá implementace předpovědi, asi 15 min na dvou jádrech) a jeho log; `provenance/` obsahuje skripty tak, jak běžely (viz `code/provenance/README.md`) |

## Reprodukce

```
pip install -r requirements.txt
python3 -I code/independent_check_v.py .        # přepočítá hlavní čísla z hodnot ve středu
python3 -I code/figures_v.py data/results/study_v_results.json figures
python3 -I code/weighted_density_independent.py . 2  # asi 15 min; opravená předpověď od začátku
```

Samotné hodnoty a nuly spočetl PARI/GP 2.17.2 (cypari2) skripty z `code/provenance/values_and_zeros/` v macOS s 9 procesy (nuly D100b, D100c, S7 a S9L 7,0, 8,2, 6,1 a 19,7 hodiny); nuly L-funkcí přiložené nejsou, manifesty se SHA-256 každého bloku ano. Ostatní skripty v `code/provenance/` obsahují cesty počítače, na kterém běžely, a jsou přiložené jako záznam.

Text, obrázky a data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`); kód: MIT (`LICENSE`). České shrnutí licencí: `LICENCE.cs.md`.
