# AntiMode — Submission Tracker
**Target**: Journal of Nonparametric Statistics (GNST, Taylor & Francis)
**Manuscript**: `submission/antimode_GNST_manuscript.tex`
**Author**: Samir Orujov — sorujov@ada.edu.az — ADA University / ICTA / UNEC

---

## Pipeline

| Command | Effect |
|---|---|
| `make fast` | R=50 quick run → tables → PDF (~15 min) |
| `make simulate` | R=500 publication run (~several hours) |
| `make tables` | CSV → .tex only (after simulate) |
| `make pdf` | pdflatex ×2 only |
| `make clean` | remove aux files + result CSVs |

---

## Completed

- [x] Manuscript restructured to `interact.cls` (Taylor & Francis GNST template)
- [x] Abstract trimmed to ≤150 words; keywords reduced to 5; MSC codes added
- [x] **Asymptotic variance V** corrected — was divergent (contained n,h); now `V = ‖K'‖²q(p*)⁵/(c₀³γ⁴)`, a pure constant
- [x] **Bias Bn** corrected — now expressed via q'''(p*) with full chain-rule derivation `f'''(x*) = -q'''(p*)/q(p*)⁴` in Appendix A
- [x] **Theorem 2.3** renamed "Equivalent characterisations" — proves (i)↔(ii)↔(iii) only (monotone transformations of f); Fisher information characterisation moved to **Corollary 2.4** with explicit conditions (first-order condition of σ²_Q, or p*=½)
- [x] Fang & Santos (2019) citation added and in bibliography
- [x] Bootstrap proof strengthened — Hadamard differentiability conditions stated explicitly
- [x] Theorem 4.3 proof sketch added (previously absent)
- [x] `\usepackage{mathtools}` added (required for `\coloneqq`)
- [x] WolfsonChamberlain2025 → WolfsonChamberlain2026 fixed everywhere
- [x] "income (or quantile) levels" → "quantile levels" in algorithm section
- [x] Remark on practical bias elimination added (undersmoothing vs plug-in)
- [x] Table 3 restructured with explicit FDR column; n=250 footnote added
- [x] Silverman column **removed** from Table 4 (was fabricated approximation; could produce negative values)
- [x] Project reorganised: `src/`, `data/`, `results/`, `submission/`
- [x] Pipeline: `src/simulate.py` → `results/*.csv` → `src/csv_to_tex.py` → `submission/tables/*.tex` → PDF
- [x] `Makefile` created with targets: all, fast, simulate, tables, pdf, clean
- [x] Author info: Samir Orujov, sorujov@ada.edu.az, ADA University / ICTA / UNEC
- [x] Funding statement: "no specific grant from any funding agency"
- [x] S&P 500 data source: Robert Shiller's public dataset (Yale)

---

## Outstanding — Ordered by Priority

### CRITICAL (blocks submission)

- [ ] **Run full simulation on cluster** (`python src/simulate.py --full`, R=500 B=299)
  - Local run was too slow (single-threaded, ~8.4M KDE calls); moving to HPC cluster
  - After completion: copy results/*.csv back, run `make tables && make pdf`
  - Verify: RMSE slopes -0.39 to -0.42; detection >90% at n=500 for A2; coverage →95% at n=1000
  - If weak DGPs (A1, A5) show non-monotone detection, add a note in Results explaining that the fixed prominence threshold (τ=0.05) correctly rejects shallow anti-modes — by design, not a code error
  - **Consider adding `multiprocessing.Pool` to simulate.py before cluster run** to use all cores

- [ ] **Fix body text after simulation** — update claims in Results section (lines ~596–608) to match actual table values once full simulation runs

- [x] **Fix Silverman body text** — removed "Silverman's bandwidth test" from line 611 power-comparison paragraph; now correctly describes Anti-mode test vs Dip test vs Z-Dip only

### IMPORTANT (needed before submission)

- [x] **Write cover letter** — `submission/cover_letter.tex`; states topical fit, novelty (first systematic anti-mode theory, n^{2/5} CLT, bootstrap), companion paper disclosure (EcLet), no funding, no competing interests

- [x] **Data repository** — public repo at https://github.com/sorujov/Antimode-JNPS; URL filled in data availability statement and cover letter

- [x] **Companion paper cross-citation** — added to Remark 2.1 (Lorenz corollary); cites `Orujov2025EcLet`; bibliography entry added

### MODERATE (polish)

- [x] **Acknowledgements** — replaced "[Omitted for blind review]" with thanks to editor/referees + funding statement

- [x] **Word count** — updated to ~2,900 words (texcount: 2782 text + 98 headers; 357 inline + 17 displayed math not counted)

- [ ] **ZDip citation** — belongs to other authors (arXiv:2511.01705); "[Authors omitted]" must be replaced with real author names from the arXiv page. Also verify the hardcoded Z-Dip power values in `src/csv_to_tex.py` (`zdip = {250: 0.611, 500: 0.836, 1000: 0.958, 5000: 1.000}`) match what that paper reports

- [ ] **DGP B3 note** — if full simulation shows non-monotone detection for B3 (Weibull mixture), add a sentence noting the Weibull barrier is shallow and below the prominence threshold at large n; report conditional-on-detection statistics with a caveat

### MINOR (final check)

- [x] **British spelling** — grep found no American -ize variants or spelling errors; "Centers" on line 729 is a proper noun (NHANES) and unchanged
- [x] **All `\ref`, `\eqref`, `\citep` resolve** — compiled twice; zero errors, zero undefined citations
- [x] **`\url` package** — `\usepackage{hyperref}` confirmed; `\url{...}` in data availability compiles without errors

---

## Submission Checklist

- [ ] `make simulate` completed (R=500, B=299)
- [ ] `make pdf` produces clean compile (zero errors, zero undefined refs)
- [ ] Cover letter written
- [ ] Data repository live with DOI
- [ ] Author, affiliation, email confirmed (single-blind — everything visible)
- [ ] Funding statement present
- [ ] Data availability statement complete with URL
- [ ] Word count updated
- [ ] British spelling checked
- [ ] ZDip citation corrected with real author names
- [ ] Upload to Taylor & Francis portal: PDF + ZIP of LaTeX source
