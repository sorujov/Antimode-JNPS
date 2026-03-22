#!/usr/bin/env python3
"""
Monte Carlo Simulation — parallelised with multiprocessing.Pool
  "Anti-Modes, Distributional Gaps, and the Fisher Information Profile"
  Target: Journal of Nonparametric Statistics

Produces CSVs in results/:
  results/table1_bias_rmse.csv
  results/table2_rates.csv
  results/table3_coverage.csv
  results/table4_power.csv

Usage (run from project root):
  python src/simulate.py          # quick mode  (R=50,  B=49)
  python src/simulate.py --full   # publication (R=500, B=299)

After running, call `python src/csv_to_tex.py` to write LaTeX table files,
or simply run `make tables` which does both steps.

Requirements: numpy, pandas, scipy
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
from pathlib import Path
from multiprocessing import Pool
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):          # graceful fallback if tqdm missing
        return iterable

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--full",     action="store_true", help="R=500, B=299 (publication run)")
parser.add_argument("--chunk-id", type=int, default=None, help="0-based chunk index (array jobs)")
parser.add_argument("--n-chunks", type=int, default=10,  help="total number of chunks (default 10)")
parser.add_argument("--dgp",      type=str, default=None, help="run only this DGP (e.g. A2)")
ARGS = parser.parse_args()

FULL     = ARGS.full
R        = 500 if FULL else 50
B        = 299 if FULL else 49
NS       = [250, 500, 1000, 5000]
DELTA    = 0.05
TAU      = 0.05
KAPPA    = 0.05
SEED     = 42
CHUNKED  = ARGS.chunk_id is not None
CHUNK_ID = ARGS.chunk_id if CHUNKED else 0
N_CHUNKS = ARGS.n_chunks

ROOT        = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)

NCPUS = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))

# ── DGP sampler ──────────────────────────────────────────────────────
def sample(name, n, rng):
    """Returns (x, [true anti-mode percentiles])."""
    if name == "A1":
        z = rng.binomial(1, .5, n)
        return np.where(z, rng.normal(0, 1, n), rng.normal(2.6, 1, n)), [0.50]
    if name == "A2":
        z = rng.binomial(1, .5, n)
        return np.where(z, rng.normal(0, 1, n), rng.normal(2.8, 1, n)), [0.50]
    if name == "A3":
        z = rng.binomial(1, .5, n)
        return np.where(z, rng.normal(0, 1, n), rng.normal(3, 1, n)), [0.50]
    if name == "A4":
        z = rng.binomial(1, .3, n)
        return np.where(z, rng.normal(0, 1, n), rng.normal(3, 1, n)), [0.27]
    if name == "A5":
        z = rng.binomial(1, .5, n)
        return np.where(z, rng.normal(0, 1, n), rng.normal(4, 2, n)), [0.60]
    if name == "A6":
        return rng.normal(0, 1, n), []
    if name == "B1":
        z = rng.binomial(1, .5, n)
        return np.where(z, stats.skewnorm.rvs(5, 0, 1, n, random_state=rng),
                        stats.skewnorm.rvs(-5, 3, 1, n, random_state=rng)), [0.50]
    if name == "B2":
        z = rng.binomial(1, .5, n)
        return np.where(z, stats.t.rvs(3, 0, 1, n, random_state=rng),
                        stats.t.rvs(3, 2.5, 1, n, random_state=rng)), [0.50]
    if name == "B3":
        z = rng.binomial(1, .5, n)
        return np.where(z, stats.weibull_min.rvs(1.5, scale=1, size=n, random_state=rng),
                        stats.weibull_min.rvs(3.0, scale=3, size=n, random_state=rng)), [0.55]
    if name == "C1":
        z = rng.choice(3, n, p=[1/3]*3)
        return np.where(z == 0, rng.normal(-3, .8, n),
               np.where(z == 1, rng.normal(0, .8, n),
                        rng.normal(3, .8, n))), [0.33, 0.67]
    if name == "C2":
        z = rng.choice(4, n, p=[.25]*4)
        mx = [-4.5, -1.5, 1.5, 4.5]
        return np.array([rng.normal(mx[zi], .7) for zi in z]), [0.25, 0.50, 0.75]
    if name == "C3":
        u = rng.uniform(0, 1, n)
        return np.sign(u - .5) * (np.abs(u - .5) ** .3) * 4, [0.50]
    raise ValueError(f"Unknown DGP: {name}")


# ── Anti-mode estimator ───────────────────────────────────────────────
def detect(x):
    n   = len(x)
    h   = 0.6 * 1.06 * x.std() * n**(-0.2)
    pg  = np.linspace(DELTA, 1 - DELTA, 600)
    qg  = np.quantile(x, pg)
    f   = np.clip(gaussian_kde(x, bw_method=h / x.std()).evaluate(qg), 1e-10, None)
    qd  = 1.0 / f
    pk  = [i for i in range(1, len(qd) - 1)
           if qd[i] > qd[i-1] and qd[i] > qd[i+1]]
    if not pk:
        return []
    out = []
    for i in pk:
        L  = qd[:i].min()   if i > 0          else qd[i]
        RR = qd[i+1:].min() if i < len(qd)-1  else qd[i]
        valley = max(L, RR)
        pr = qd[i] - valley
        if pr >= TAU * valley:
            out.append((pg[i], pr))
    out.sort(key=lambda t: t[1], reverse=True)
    merged = []
    for p, pr in out:
        if all(abs(p - pm) > KAPPA for pm, _ in merged):
            merged.append((p, pr))
    return [p for p, _ in merged]


# ── Bootstrap CI ──────────────────────────────────────────────────────
def bci(x, ph, rng):
    bt = []
    for _ in range(B):
        xb = rng.choice(x, len(x), replace=True)
        a  = detect(xb)
        if a:
            bt.append(min(a, key=lambda p: abs(p - ph)))
    if len(bt) < 5:
        return np.nan, np.nan
    return np.percentile(bt, 2.5), np.percentile(bt, 97.5)


# ── Unimodality test ──────────────────────────────────────────────────
def am_test(x, rng):
    n   = len(x)
    pg  = np.linspace(DELTA, 1 - DELTA, 300)
    h   = 0.6 * 1.06 * x.std() * n**(-0.2)
    qg  = np.quantile(x, pg)
    f   = np.clip(gaussian_kde(x, bw_method=h / x.std()).evaluate(qg), 1e-10, None)
    qd  = 1.0 / f
    mu0, s0 = x.mean(), x.std()
    qn  = s0 / stats.norm.pdf(stats.norm.ppf(pg))
    T   = np.max(np.clip(qd - qn, 0, None))
    null = []
    for _ in range(B):
        xb      = rng.normal(mu0, s0, n)
        mb, sb  = xb.mean(), xb.std()
        hb      = 0.6 * 1.06 * sb * n**(-0.2)
        fb      = np.clip(gaussian_kde(xb, bw_method=hb / sb).evaluate(np.quantile(xb, pg)), 1e-10, None)
        qdb     = 1.0 / fb
        qnb     = sb / stats.norm.pdf(stats.norm.ppf(pg))
        null.append(np.max(np.clip(qdb - qnb, 0, None)))
    return T > np.percentile(null, 95)


def dip_test(x):
    import diptest
    _, pval = diptest.diptest(x)
    return pval < 0.05


# Null-distribution parameters for the dip statistic under Uniform[0,1].
# Estimated via 5,000 Monte Carlo replications (seed 42) for each N.
# Used to compute the Z-Dip statistic of Di Martino et al. (2025).
_ZDIP_MU = {250: 0.02364148, 500: 0.01698986, 1000: 0.01214315, 5000: 0.00546081}
_ZDIP_SD = {250: 0.00506343, 500: 0.00362469, 1000: 0.00260220, 5000: 0.00114250}

def zdip_test(x):
    """Z-Dip test (Di Martino et al. 2025) at nominal 5% level.
    Standardises the dip statistic against its null distribution under
    Uniform[0,1] and rejects when Z-Dip > 1.975 (universal threshold)."""
    import diptest
    n = len(x)
    dip_stat, _ = diptest.diptest(x)
    mu = _ZDIP_MU.get(n)
    sd = _ZDIP_SD.get(n)
    if mu is None or sd is None or sd == 0:
        # Fallback for N not in lookup: use dip p-value at 5%
        _, pval = diptest.diptest(x)
        return pval < 0.05
    z = (dip_stat - mu) / sd
    return z > 1.975


# ── Single replication (worker function) ─────────────────────────────
def _one_rep(args):
    dgp, n, pt, kt, rep_seed = args
    rng  = np.random.default_rng(rep_seed)
    x, _ = sample(dgp, n, rng)
    am   = detect(x)
    kh   = len(am)

    if dgp == "A6":
        return dict(bs=np.nan, sq=np.nan, cv=np.nan,
                    dt=1 if kh == 0 else 0,
                    fd=1 if kh > 0  else 0,
                    tam=1   if am_test(x, rng)  else 0,
                    tdip=1  if dip_test(x)      else 0,
                    tzdip=1 if zdip_test(x)     else 0)
    else:
        dt   = 1 if kh == kt else 0
        bs_v = sq_v = cv_v = np.nan
        if kh > 0 and kt > 0:
            ph   = min(am, key=lambda p: abs(p - pt[0]))
            e    = ph - pt[0]
            bs_v = e
            sq_v = e**2
            lo, hi = bci(x, ph, rng)
            if not np.isnan(lo):
                cv_v = 1.0 if lo <= pt[0] <= hi else 0.0
        return dict(bs=bs_v, sq=sq_v, cv=cv_v,
                    dt=dt, fd=0,
                    tam=1   if am_test(x, rng)  else 0,
                    tdip=1  if dip_test(x)      else 0,
                    tzdip=1 if zdip_test(x)     else 0)


# ── Main ──────────────────────────────────────────────────────────────
DGPS_ALL = ["A6", "A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "C1", "C2", "C3"]
DGPS = [ARGS.dgp] if ARGS.dgp else DGPS_ALL


def _aggregate(reps, dgp):
    """Aggregate a list of per-replication dicts into summary stats."""
    bs    = [r["bs"]    for r in reps]
    sq    = [r["sq"]    for r in reps]
    cv    = [r["cv"]    for r in reps if not np.isnan(r["cv"])]
    dt    = [r["dt"]    for r in reps]
    fd    = [r["fd"]    for r in reps]
    tam   = [r["tam"]   for r in reps]
    tdip  = [r["tdip"]  for r in reps]
    tzdip = [r["tzdip"] for r in reps]
    return dict(
        Bias       = round(np.nanmean(bs),  4) if dgp != "A6" else None,
        RMSE       = round(float(np.sqrt(np.nanmean(sq))), 4) if dgp != "A6" else None,
        Coverage   = round(np.mean(cv) * 100, 1) if cv else None,
        Detection  = round(np.mean(dt) * 100, 1),
        FDR        = round(np.mean(fd) * 100, 1),
        Power_AM   = round(np.mean(tam)   * 100, 1) if dgp != "A6" else None,
        Power_Dip  = round(np.mean(tdip)  * 100, 1) if dgp != "A6" else None,
        Power_ZDip = round(np.mean(tzdip) * 100, 1) if dgp != "A6" else None,
    )


def main():
    chunk_tag = f"_chunk{CHUNK_ID:02d}" if CHUNKED else ""
    mode_str  = f"chunk {CHUNK_ID+1}/{N_CHUNKS}" if CHUNKED else "full"
    print(f"R={R}, B={B}, mode={mode_str}, workers={NCPUS}")

    # Replication index range for this chunk
    reps_per_chunk = max(1, R // N_CHUNKS)
    r_start = CHUNK_ID * reps_per_chunk
    r_end   = R if (CHUNK_ID == N_CHUNKS - 1) else r_start + reps_per_chunk
    r_range = range(r_start, r_end)

    tab1, tab3, tab4 = [], [], []
    raw_rows = []    # saved in chunk mode for later merging

    cells = [(dgp, n) for dgp in DGPS for n in NS]

    with Pool(processes=NCPUS) as pool:
        for dgp, n in tqdm(cells, desc="cells", unit="cell"):
            _, pt = sample(dgp, 1000, np.random.default_rng(SEED))
            kt    = len(pt)

            args = [(dgp, n, pt, kt, SEED + r * 997 + NS.index(n) * 13)
                    for r in r_range]

            reps = list(tqdm(
                pool.imap(_one_rep, args),
                total=len(args), desc=f"{dgp} n={n}", leave=False
            ))

            agg = _aggregate(reps, dgp)
            tab1.append({"DGP": dgp, "n": n, "Bias": agg["Bias"], "RMSE": agg["RMSE"]})
            tab3.append({"DGP": dgp, "n": n, "Coverage": agg["Coverage"],
                         "Detection": agg["Detection"], "FDR": agg["FDR"]})
            tab4.append({"DGP": dgp, "n": n, "Power_AM": agg["Power_AM"],
                         "Power_Dip": agg["Power_Dip"],
                         "Power_ZDip": agg["Power_ZDip"]})
            print(f"  {dgp} n={n}: bias={agg['Bias']} rmse={agg['RMSE']} "
                  f"det={agg['Detection']}% cov={agg['Coverage']}%")

            if CHUNKED:
                for rep_idx, rep in zip(r_range, reps):
                    raw_rows.append({"chunk": CHUNK_ID, "rep": rep_idx,
                                     "DGP": dgp, "n": n, **rep})

    # In chunk mode: save raw replication data for merging, plus partial aggregates
    if CHUNKED:
        pd.DataFrame(raw_rows).to_csv(
            RESULTS_DIR / f"raw_reps{chunk_tag}.csv", index=False)
        print(f"Raw reps saved → results/raw_reps{chunk_tag}.csv")

    # Table 2: convergence rates with bootstrap CIs
    tab2    = []
    log_ns  = np.log(NS)
    for dgp in [d for d in DGPS if d != "A6"]:
        rmses = [next((r["RMSE"] for r in tab1 if r["DGP"] == dgp and r["n"] == n), None)
                 for n in NS]
        valid = [(ln, np.log(r)) for ln, r in zip(log_ns, rmses) if r and r > 0]
        if len(valid) < 2:
            continue
        xs    = [v[0] for v in valid]
        ys    = [v[1] for v in valid]
        slope = float(np.polyfit(xs, ys, 1)[0])
        rng_b    = np.random.default_rng(SEED)
        slopes_b = []
        for _ in range(2000):
            idx  = rng_b.choice(len(valid), len(valid), replace=True)
            xs_b = [xs[i] for i in idx]
            ys_b = [ys[i] for i in idx]
            if len(set(xs_b)) < 2:
                continue
            slopes_b.append(float(np.polyfit(xs_b, ys_b, 1)[0]))
        lo_ci = round(float(np.percentile(slopes_b, 2.5)),  3) if slopes_b else None
        hi_ci = round(float(np.percentile(slopes_b, 97.5)), 3) if slopes_b else None
        tab2.append({"DGP": dgp, "Slope": round(slope, 3),
                     "CI_lo": lo_ci, "CI_hi": hi_ci})

    pd.DataFrame(tab1).to_csv(RESULTS_DIR / f"table1_bias_rmse{chunk_tag}.csv", index=False)
    pd.DataFrame(tab2).to_csv(RESULTS_DIR / f"table2_rates{chunk_tag}.csv",     index=False)
    pd.DataFrame(tab3).to_csv(RESULTS_DIR / f"table3_coverage{chunk_tag}.csv",  index=False)
    pd.DataFrame(tab4).to_csv(RESULTS_DIR / f"table4_power{chunk_tag}.csv",     index=False)
    print(f"\nResults saved to {RESULTS_DIR}/ (tag='{chunk_tag}')")
    print(pd.DataFrame(tab2).to_string(index=False))


if __name__ == "__main__":
    main()
