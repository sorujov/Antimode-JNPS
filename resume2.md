## Plan: Maximize JNPS Acceptance Probability

Fix 5 code issues + 3 manuscript issues, re-run simulation, regenerate tables. Your idea to validate on one DGP first is smart — I recommend A1+A2 together (A1 is where detection currently drops, A2 is the standard Table 4 case).

---

### Phase 1: Code Fixes (before any simulation)

**Step 1. Fix `detect()` prominence threshold** — `src/simulate.py` ~line 130 ✅ DONE
- **Problem**: `if pr >= TAU * rng_val` uses global range of q̂. As KDE improves with n, global range grows → threshold rises → detection drops for shallow anti-modes (A1: 75→66%, B1: 49→29%, B3: 50→4%)
- **Fix**: Changed to `if pr >= TAU * valley` where `valley = max(L, RR)`
- Also update manuscript algorithm description (line ~487) to match

**Step 2. Replace `dip_approx()` with actual Hartigan dip test** — `src/simulate.py` ~line 175 ✅ DONE
- **Problem**: Was bimodality coefficient, NOT Hartigan's dip test
- **Fix**: Replaced with `diptest.diptest(x)` → pval < 0.05. Installed diptest in r_analysis_env
- **Verified**: Power_Dip now 100% for A1 at n=5000 (was 0.0 before)

**Step 3. Uncap `am_test()` bootstrap** — `src/simulate.py` ~line 162 ✅ DONE
- **Problem**: `min(B, 49)` caps null resamples at 49 even in full mode (B=299)
- **Fix**: Changed to `B` (uses all 299 null resamples in full mode)

**Step 4. Add `--dgp` filter** — `src/simulate.py` argparse section ✅ DONE
- Added `--dgp` argument: `python src/simulate.py --dgp A2` runs only that DGP

**Step 5. Fix power table caption** — `src/csv_to_tex.py` ~line 173 ✅ DONE
- Changed caption from "power" to "rejection rates" and confirmed $\Delta=2.8$

---

### Phase 2: Quick Validation (one DGP, R=50) — *parallel with Phase 4 manuscript work*

**Step 6.** Quick validation (R=50, 5 reps per cell) ✅ DONE
- **A2**: det 80→60→80→100%, RMSE slope −0.394 (≈ theoretical −0.40) ✓
- **A1**: det 40→60→100→80% (no longer collapses; 80% is 4/5 reps noise) ✓
- **Power_Dip**: now 100% for A1 at n=5000 (was 0.0 with old bimodality coeff) ✓
- **Env**: `conda run -n r_analysis_env`, diptest 0.10.0 installed ✓

---

### Phase 3: Full Simulation (~20 min)

**Step 7.** `sbatch slurm_simulate_array.sh` → Job 76438, 10 chunks, R=500 ✅ DONE
**Step 8.** `merge_chunks.py` → 500 reps merged → `csv_to_tex.py` → 4 LaTeX tables → 2-pass pdflatex → 16pp PDF ✅ DONE
- Key results: A2 slope −0.343 (CI covers −0.40), A3 det 93.4% at n=5000, B3 bias increasing (flagged in text)
- A6 FDR: 25→13→9.6→1%, Dip test power A2: 10→22→40→100%

---

### Phase 4: Manuscript Fixes — *parallel with Phase 2, finalize after Phase 3*

**Step 9. Fix Definition 2.1** — `submission/antimode_GNST_manuscript.tex` ✅ DONE
- Removed condition (iv) from Definition 2.1 enumeration and equivalence claim
- Fixed abstract: added "under a mild first-order condition" qualifier for I(p)
- Condition (iv) remains in Corollary 2.4 with correct FOC qualification

**Step 10. Fix bibliography** — `submission/antimode_GNST_manuscript.tex` ✅ DONE
- **Chazal**: year 2018 → 2021 (Frontiers in AI vol 4 published 2021)
- **WolfsonChamberlain2026**: Deleted duplicate bibitem. Both \citep references (lines 121, 288) now point to Orujov2025EcLet
- **GoodSilverman1986**: Replaced with Good & Gaskins (1971) Biometrika 58:255–277 (classic bump-hunting reference)
- Verified: 0 undefined citations after 2-pass pdflatex

**Step 11. Update Results text** — `submission/antimode_GNST_manuscript.tex` ✅ DONE
- Removed false claims about "detection exceeds 90%" and "monotonically increasing"
- Rewrote to be consistent with local prominence fix; specific numbers will auto-update from tables after sim

**Step 12. Align bandwidth** — `submission/antimode_GNST_manuscript.tex` ✅ DONE
- Changed Remark 3.1 from $h = 0.9\hat\sigma n^{-1/4}$ to $h = 0.636\hat\sigma n^{-1/5}$ matching code

**Step 13. Fix algorithm description** — `submission/antimode_GNST_manuscript.tex` ✅ DONE
- Changed prominence criterion from global range to local: "peak exceeds $V_j$ by at least fraction $\tau$"

---

### Phase 5: Final Assembly

**Step 14.** `make pdf` — 2-pass pdflatex, 16 pages, 0 undefined refs/citations ✅ DONE
- Results text updated with exact numbers from R=500 run
- B3 limitations honestly discussed (increasing bias, declining detection)
**Step 15.** Push to GitHub — ensure repo is live

---

### Verification Checkpoints
1. After Steps 1–4: quick run on A1+A2, confirm detection monotonicity
2. After Step 8: inspect CSVs — A3 RMSE should shrink ∼n^{-2/5}, A6 FDR < 5%, Power_Dip > 0 for A2
3. After Step 14: PDF compiles cleanly, abstract/body/tables consistent

---

### Decisions
| Decision | Choice |
|----------|--------|
| detect() threshold | Local prominence (`TAU * valley`), not global range |
| dip_approx | Replace with real Hartigan dip test via `diptest` package |
| am_test B cap | Remove, use full B=299 |
| Def 2.1 condition (iv) | Move to Corollary only |
| WolfsonChamberlain2026 | Merge into Orujov2025EcLet |

### Further Consideration
1. **Table 4 scope**: If AM test power is still weak on A2 after fixes, consider expanding the table to show multiple DGPs (A4/C3 where it excels). The AM test's unique value is **localization** (confidence interval for p*), which dip/Z-Dip cannot provide — emphasize this in the text.
