# Anti-Modes, Distributional Gaps, and the Fisher Information Profile

**Nonparametric Estimation and Inference**

Submitted to the *Journal of Nonparametric Statistics* (Taylor & Francis)

---

## What is an anti-mode?

A **mode** is where a density peaks. An **anti-mode** is where it dips — the sparse zone between two population clusters. This paper provides the first complete nonparametric theory for locating, estimating, and testing anti-modes.

<p align="center">
<em>Where the density is thinnest, estimation uncertainty is highest, and the Fisher information is lowest.</em>
</p>

## Key idea

We introduce the **Fisher information profile**

$$I(p) = \frac{f(Q(p))^2}{p(1-p)} = \frac{1}{\sigma^2_Q(p)}$$

and show that anti-modes are simultaneously:

| Characterisation | Condition at $p^*$ |
|---|---|
| Density | local minimum of $f$ |
| Quantile density | local maximum of $q(p) = 1/f(Q(p))$ |
| Estimation variance | local maximum of $\sigma^2_Q(p)$ |
| Fisher information | local minimum of $I(p)$ |

These four views are equivalent — and the last is new.

## Results

| Result | Rate / Guarantee |
|---|---|
| Estimator consistency | $n^{2/5}$ |
| Asymptotic normality | $n^{2/5}(\hat{p}^* - p^*) \xrightarrow{d} N(0, V)$ |
| Uniform convergence of $\hat{I}(p)$ | $O_p((nh)^{-1/2} + h^2)$ |
| Bootstrap CI coverage | $\geq 95\%$ at moderate $n$ |
| AMSE-optimal rate | $n^{-4/5}$ (minimax) |

## Empirical illustration

Applied to the **Old Faithful geyser** eruption data ($n = 272$):

- Anti-mode detected at $\hat{p}^* = 0.357$ → 3.0 minutes (gap between short and long eruptions)
- 95% bootstrap CI: $[2.30,\ 3.72]$ minutes
- All tests reject unimodality ($p < 0.0001$)

## Simulation

12 DGPs across three families (Gaussian mixtures, heavy-tailed/skew mixtures, multi-modal), evaluated at $n \in \{250, 500, 1000, 5000\}$ with $R = 500$ replications and $B = 299$ bootstrap draws.

## Repository structure

```
├── submission/
│   ├── antimode_GNST_manuscript.tex   # Manuscript (LaTeX source)
│   ├── cover_letter.tex               # Cover letter
│   └── tables/                        # Auto-generated LaTeX tables
├── src/
│   ├── simulate.py                    # Monte Carlo simulation engine
│   ├── empirical_faithful.py          # Old Faithful empirical illustration
│   ├── csv_to_tex.py                  # Results → LaTeX table converter
│   └── merge_chunks.py               # SLURM array job merger
├── slurm_simulate.sh                  # SLURM batch script (full run)
├── slurm_simulate_array.sh            # SLURM array job script
└── Makefile                           # make simulate / make tables
```

## Reproduce

```bash
# Quick run (R=50, ~10 min on 32 cores)
python src/simulate.py

# Full publication run (R=500, submit to SLURM)
sbatch slurm_simulate.sh

# Generate LaTeX tables from results
python src/csv_to_tex.py
```

## Citation

```bibtex
@article{Orujov2026antimode,
  title   = {Anti-Modes, Distributional Gaps, and the Fisher Information
             Profile: Nonparametric Estimation and Inference},
  author  = {Orujov, Samir},
  journal = {Journal of Nonparametric Statistics},
  year    = {2026},
  note    = {Submitted}
}
```

## Author

**Samir Orujov**
ADA University / ICTA / UNEC, Baku, Azerbaijan
Université Bretagne Sud, Vannes, France
[sorujov@ada.edu.az](mailto:sorujov@ada.edu.az) · [ORCID](https://orcid.org/0009-0004-9708-2109)
