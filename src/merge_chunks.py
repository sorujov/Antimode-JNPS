#!/usr/bin/env python3
"""
Merge raw replication CSVs produced by array-job chunks into final result tables.

Usage (run from project root after all array chunks complete):
  python src/merge_chunks.py
  make merge-chunks        # equivalent

Reads:   results/raw_reps_chunk*.csv
Writes:  results/table{1-4}_*.csv  (same format as simulate.py non-chunked output)
Then:    run `make tables` to regenerate LaTeX tables from the merged CSVs.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from glob import glob

ROOT        = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"

NS   = [250, 500, 1000, 5000]
DGPS = ["A6", "A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "C1", "C2", "C3"]
SEED = 42


def load_all_raw():
    files = sorted(glob(str(RESULTS_DIR / "raw_reps_chunk*.csv")))
    if not files:
        raise SystemExit("No raw_reps_chunk*.csv files found in results/. "
                         "Run the array job first.")
    print(f"Found {len(files)} chunk files: {[Path(f).name for f in files]}")
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def aggregate(df, dgp, n):
    sub  = df[(df.DGP == dgp) & (df.n == n)]
    bs   = sub["bs"].values
    sq   = sub["sq"].values
    cv   = sub["cv"].dropna().values
    dt   = sub["dt"].values
    fd   = sub["fd"].values
    tam   = sub["tam"].values
    tdip  = sub["tdip"].values
    tzdip = sub["tzdip"].values
    return dict(
        Bias       = round(float(np.nanmean(bs)),  4) if dgp != "A6" else None,
        RMSE       = round(float(np.sqrt(np.nanmean(sq))), 4) if dgp != "A6" else None,
        Coverage   = round(float(np.mean(cv)) * 100, 1) if len(cv) > 0 else None,
        Detection  = round(float(np.mean(dt)) * 100, 1),
        FDR        = round(float(np.mean(fd)) * 100, 1),
        Power_AM   = round(float(np.mean(tam))   * 100, 1) if dgp != "A6" else None,
        Power_Dip  = round(float(np.mean(tdip))  * 100, 1) if dgp != "A6" else None,
        Power_ZDip = round(float(np.mean(tzdip)) * 100, 1) if dgp != "A6" else None,
    )


def build_tab2(tab1):
    log_ns   = np.log(NS)
    tab2     = []
    rng_b    = np.random.default_rng(SEED)
    for dgp in [d for d in DGPS if d != "A6"]:
        rmses = [next((r["RMSE"] for r in tab1 if r["DGP"] == dgp and r["n"] == n), None)
                 for n in NS]
        valid = [(ln, np.log(r)) for ln, r in zip(log_ns, rmses) if r and r > 0]
        if len(valid) < 2:
            continue
        xs    = [v[0] for v in valid]
        ys    = [v[1] for v in valid]
        slope = float(np.polyfit(xs, ys, 1)[0])
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
    return tab2


def main():
    df = load_all_raw()
    total_reps = df["rep"].nunique()
    print(f"Total unique replications loaded: {total_reps}")

    tab1, tab3, tab4 = [], [], []
    for dgp in DGPS:
        for n in NS:
            agg = aggregate(df, dgp, n)
            tab1.append({"DGP": dgp, "n": n, "Bias": agg["Bias"], "RMSE": agg["RMSE"]})
            tab3.append({"DGP": dgp, "n": n, "Coverage": agg["Coverage"],
                         "Detection": agg["Detection"], "FDR": agg["FDR"]})
            tab4.append({"DGP": dgp, "n": n, "Power_AM": agg["Power_AM"],
                         "Power_Dip": agg["Power_Dip"],
                         "Power_ZDip": agg["Power_ZDip"]})

    tab2 = build_tab2(tab1)

    pd.DataFrame(tab1).to_csv(RESULTS_DIR / "table1_bias_rmse.csv", index=False)
    pd.DataFrame(tab2).to_csv(RESULTS_DIR / "table2_rates.csv",     index=False)
    pd.DataFrame(tab3).to_csv(RESULTS_DIR / "table3_coverage.csv",  index=False)
    pd.DataFrame(tab4).to_csv(RESULTS_DIR / "table4_power.csv",     index=False)
    print("Merged tables saved to results/")
    print(pd.DataFrame(tab2).to_string(index=False))
    print("\nNext: run `make tables` to regenerate LaTeX files.")


if __name__ == "__main__":
    main()
