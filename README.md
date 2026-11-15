# Central values of quadratic Dirichlet L-functions for every fundamental discriminant up to 10⁷ and at 10⁸, 10⁹ and 10¹⁰ (study V)

**David Hess (MAYAN ALFA V1008 LAB v.34.85b) · 2026** · [Česká verze](README.cs.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24256378.svg)](https://doi.org/10.5281/zenodo.24256378) Zenodo: [doi:10.5281/zenodo.24256378](https://doi.org/10.5281/zenodo.24256378)

Fifth study of the series, after *Berry saturation and the statistics of the Riemann zeros* (study I, [doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328)), *Beyond the diagonal approximation* (study II, [doi:10.5281/zenodo.23928017](https://doi.org/10.5281/zenodo.23928017)), *Low-lying zeros of 795,381 quadratic Dirichlet L-functions* (study III, [doi:10.5281/zenodo.23991165](https://doi.org/10.5281/zenodo.23991165)) and *Three-point correlation of the Riemann zeros at finite height* (study IV, [doi:10.5281/zenodo.24094318](https://doi.org/10.5281/zenodo.24094318)). It measures the values L(½, χ_D) of quadratic Dirichlet L-functions at the central point and connects them with the low-lying zeros of the same L-functions (study III). Data: every fundamental discriminant 3 ≤ |D| ≤ 10⁷ (6,079,285 values, PARI/GP), all 1,215,839 discriminants with 10⁸ < |D| ≤ 10⁸ + 2·10⁶ (PARI/GP), all 486,344 with 10⁹ < |D| ≤ 10⁹ + 8·10⁵ and all 121,557 with 10¹⁰ < |D| ≤ 10¹⁰ + 2·10⁵ (own smoothed functional equation, checked against PARI on 2,000 random discriminants at 10⁹ and on all 12,153 of the first band at 10¹⁰); for 4,042,611 of them also the low zeros (the zeros of 3,161,211 L-functions were computed for this study). Most statements were tested against rules recorded before the computation (`data/rules/`).

| result | value |
|---|---|
| moments k = 1–5 against the full CFKRS polynomials, 7 bands of \|D\| ≤ 10⁷ × 2 signs × 5 (70 combinations) | largest deviation 0.65 standard errors; totals 6·10⁻⁷ – 3·10⁻⁴ relative |
| the same at 10⁸ (D100b ∪ D100c, 1,155,037 values) | largest deviation 0.57 errors |
| descriptively (no rule recorded) the same at 10⁹ (S9L, 486,344 values) | largest deviation 1.15 errors |
| leading term alone (Keating–Snaith with the arithmetic factor) | off by 9 % (k = 1) to a factor 6·10⁴ (k = 5) |
| our CFKRS polynomials against the predicted values of Alderson and Rubinstein (k = 1–5, both signs, to X = 5·10¹⁰) | 6.5·10⁻⁸ |
| difference of the first moments, odd − even: data / A₀π/2 | 1.10649 ± 0.00060 / 1.10654 |
| sums over blocks of 20,000 consecutive \|D\| against random samples of the same size (k = 1–4) | scatter 0.11–0.33 of that of independent values (random samples 0.94–1.04); falls with the width w; scaling with w/√\|D\| and with ln w/ln √\|D\| rejected |
| variance of ln L against a hybrid (arithmetic factor × USp(2N), N from the first moment, no free parameter) | within 0.26 % |
| exponent of small values, P(L < ε) ∝ ε^β | β = 1.518 ± 0.015 (USp: 3/2); 11–17 % fewer small values than the hybrid |
| smallest value; values ≤ 0 | 6.04·10⁻⁵ at D = −115,147; none |
| correlation of ln L(½) with the first zero: data / USp(2N) with N from the first moment | 0.83–0.89 / 0.81–0.89 (data above in 7 of 8 groups, up to 17 errors) |
| N_eff from the ratios conjecture: mean first zero / its shift under weighting by L(½) | explained / still 9–16 % above USp(2N_eff) |
| weighted one-level density (ratios conjecture with weight L(½), Fazzari), both signs, 10⁴–10¹⁰ | χ² 0.2–4.3 for 5 bins in 20 groups (smallest p 0.51; the errors are cautious, Section 4.6) |
| a family computed after the correction under a rule recorded before its data (S7: every 7.5·10⁶ < \|D\| ≤ 10⁷, 1,519,830 L, 6,592,004 zeros) | χ² 0.21 / 0.64 (p 0.999 / 0.986); even amplitude 0.015 ± 0.040 (the erroneous implementation would give 0.62, as its scaling law had predicted: 0.608 ± 0.045) |
| the family at 10⁹, computed after the correction under a rule recorded before its data (S9L: every 10⁹ < \|D\| ≤ 10⁹ + 8·10⁵, 486,344 L, 2,099,543 zeros; values by the own functional equation) | χ² 1.75 / 1.31 (p 0.883 / 0.934); even amplitude −0.088 ± 0.084 (the erroneous implementation would give 0.13) |
| correction before publication: our first implementation of this prediction (one wrong Γ factor) produced an apparent term of even characters falling like \|D\|^(−1/4); an independent implementation recorded beforehand found the error | difference up to 13.3 errors before, 0.005 errors after the correction; even amplitudes now 0.065 ± 0.059 (C100), −0.038 ± 0.050 (10⁸) |
| explicit formula on C100 (607,941 L) and S7 (1,519,830 L): zero side against prime side from the measured L(½), weighted − unweighted | agreement to 2.2 / 1.6 errors on C100 and 2.0 / 0.1 on S7 (even / odd) |
| L(½) at 10¹⁰: own functional equation against PARI, 12,153 values | largest relative difference 4.3·10⁻¹¹ |
| appendix: Ford–Zaharescu on 2,001,052 zeros of ζ, {αγ} mod 1 against Landau's formula | smallest p 0.66 (α = ln p/2π), uniform excluded for p = 2; controls uniform |
| independent check (`code/independent_check_v.py`, no shared code) | 12 of 12 checks passed (`code/independent_check_v.log`) |

The CFKRS polynomials, the exponent 3/2, the hybrid models and the weighted density of the ratios conjecture are known theory (conjectural except where proved); the contribution is the joint data of values and zeros, the effect of consecutive families on the moments, the comparisons without free parameters, and a test of the weighted density of the ratios conjecture for both signs on 4.0 million L-functions, including two families (at 10⁷ and at 10⁹) computed after the correction under rules recorded before their data, with a check by the explicit formula on two families. An apparent lower-order term of even characters, found with our first implementation of the prediction, was an error of that implementation; it is described in Section 4.6 of the article. Agreement is numerical support for a conjecture, not a proof; nothing here bears on the generalised Riemann hypothesis.

## Contents

| folder | what |
|---|---|
| `paper/` | the article `central_values_quadratic_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); the Czech version in `paper/cs/` |
| `figures/` | figures 1–9 in English (`_en`) and Czech (`_cz`), PNG 300 dpi and PDF |
| `data/central_values/` | L(½, χ_D) for every discriminant of Table 1 of the paper (NumPy `.npz`, arrays `D` and `L`, sorted by \|D\|; parts below 20 MB), `index.json` with counts and SHA-256 |
| `data/results/` | `study_v_results.json`: every number of the paper (moments by band, blocks, ln L, small values, first zero, hybrid products, weighted density with the corrected prediction, the families S7 and S9L, the explicit formula on C100 and S7, the earlier numbers of the erroneous implementation for the record, Ford–Zaharescu, verdicts) |
| `data/cfkrs/` | coefficients of the CFKRS polynomials Q_k(ln\|D\|), k = 1–5, both signs; the tables of Alderson and Rubinstein used for the comparison |
| `data/manifests/` | manifests of the runs (SHA-256 of every block of values and zeros, environment, script checksum) |
| `data/rules/` | the rules recorded before the computation, with dates (Czech, as written) |
| `code/` | `figures_v.py` (figures from the results file), `independent_check_v.py` and its log, `weighted_density_independent.py` (the independent implementation of the prediction, about 15 min on two cores) and its log; `provenance/` holds the scripts as they were run (see `code/provenance/README.md`) |

## Reproducing

```
pip install -r requirements.txt
python3 -I code/independent_check_v.py .        # recomputes the main numbers from the central values
python3 -I code/figures_v.py data/results/study_v_results.json figures
python3 -I code/weighted_density_independent.py . 2  # about 15 min; the corrected prediction from scratch
```

The central values and zeros themselves were computed with PARI/GP 2.17.2 (cypari2) by the scripts in `code/provenance/values_and_zeros/` on macOS with 9 processes (zeros of D100b, D100c, S7 and S9L 7.0, 8.2, 6.1 and 19.7 hours); the zeros of the L-functions are not included, but the manifests with the SHA-256 of every block are. The other scripts in `code/provenance/` contain the paths of the computer on which they ran and are given as a record.

Text, figures and data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`); code: MIT (`LICENSE`). Czech summary of the licences: `LICENCE.cs.md`.
