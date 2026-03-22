#!/usr/bin/env python3
"""
Convert simulation CSVs to LaTeX table files.

Reads:   results/table{1-4}_*.csv
Writes:  submission/tables/tab_{bias_rmse,rates,coverage,power}.tex

Usage (run from project root):
  python src/csv_to_tex.py

The LaTeX files are input{}-ed by submission/antimode_GNST_manuscript.tex,
so recompiling the PDF after running this script picks up new results.
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT       = Path(__file__).resolve().parent.parent
RESULTS    = ROOT / "results"
TABLES_OUT = ROOT / "submission" / "tables"
TABLES_OUT.mkdir(parents=True, exist_ok=True)

NS  = [250, 500, 1000, 5000]
DGPS_ALL = ["A6", "A1", "A2", "A3", "A4", "A5",
            "B1", "B2", "B3", "C1", "C2", "C3"]
NULL_DGP = "A6"   # unimodal, no true anti-mode

# Human-readable labels for each DGP
LABELS = {
    "A1": r"A1 (weak, $\Delta=2.6$)",
    "A2": r"A2 (moderate, $\Delta=2.8$)",
    "A3": r"A3 (sharp, $\Delta=3$)",
    "A4": r"A4 (unequal $\pi$, $\Delta=3$)",
    "A5": r"A5 (unequal $\sigma$, $\Delta\approx2.7$)",
    "A6": r"A6 (unimodal, size)",
    "B1": r"B1 (skew-normal)",
    "B2": r"B2 (heavy tail $t_3$)",
    "B3": r"B3 (Weibull mixture)",
    "C1": r"C1 (two anti-modes)",
    "C2": r"C2 (three anti-modes)",
    "C3": r"C3 (bathtub)",
}


def _fmt(v, decimals=3, fallback="---"):
    """Format a float; return fallback if None or NaN."""
    if v is None:
        return fallback
    try:
        if np.isnan(float(v)):
            return fallback
    except (TypeError, ValueError):
        return fallback
    return f"{float(v):.{decimals}f}"


# ── Table 1: Bias and RMSE ────────────────────────────────────────────────────
def write_tab_bias_rmse(df):
    lines = [
        r"\clearpage",
        r"\begin{table}[p]",
        r"\centering",
        r"\caption{Bias and RMSE of $\hat p^*$ (percentile-point units)",
        r"across DGPs and sample sizes; $R=500$, $B=299$.}",
        r"\label{tab:bias_rmse}",
        r"\begin{tabular}{lcccccccc}",
        r"\toprule",
        r" & \multicolumn{2}{c}{$n=250$} & \multicolumn{2}{c}{$n=500$}"
        r" & \multicolumn{2}{c}{$n=1000$} & \multicolumn{2}{c}{$n=5000$}\\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}",
        r"DGP & Bias & RMSE & Bias & RMSE & Bias & RMSE & Bias & RMSE\\",
        r"\midrule",
    ]
    for dgp in DGPS_ALL:
        row = [LABELS[dgp]]
        for n in NS:
            sub = df[(df.DGP == dgp) & (df.n == n)]
            if sub.empty or dgp == NULL_DGP:
                row += ["---", "---"]
            else:
                row.append(_fmt(sub.iloc[0]["Bias"]))
                row.append(_fmt(sub.iloc[0]["RMSE"]))
        lines.append(" & ".join(row) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TABLES_OUT / "tab_bias_rmse.tex").write_text("\n".join(lines))
    print("  Wrote tab_bias_rmse.tex")


# ── Table 2: Convergence rates ────────────────────────────────────────────────
def write_tab_rates(df):
    lines = [
        r"\clearpage",
        r"\begin{table}[p]",
        r"\centering",
        r"\caption{Empirical convergence rates: slope of",
        r"$\log(\mathrm{RMSE})$ on $\log(n)$ (theoretical: $-0.40$).}",
        r"\label{tab:rates}",
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r"DGP & Estimated slope & 95\% CI\\",
        r"\midrule",
    ]
    rate_dgps = {"A2": "A2", "A3": "A3", "B2": "B2",
                 "B3": "B3", "C1": r"C1 (barrier 1)", "C3": "C3"}
    for dgp, label in rate_dgps.items():
        sub = df[df.DGP == dgp]
        if sub.empty:
            lines.append(f"{label} & --- & ---\\\\")
        else:
            sl = _fmt(sub.iloc[0]["Slope"])
            lo = _fmt(sub.iloc[0]["CI_lo"])
            hi = _fmt(sub.iloc[0]["CI_hi"])
            lines.append(f"{label} & ${sl}$ & $[{lo},\\;{hi}]$\\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TABLES_OUT / "tab_rates.tex").write_text("\n".join(lines))
    print("  Wrote tab_rates.tex")


# ── Table 3: Coverage, detection, FDR ────────────────────────────────────────
def write_tab_coverage(df):
    ns_cov = [500, 1000, 5000]
    lines = [
        r"\clearpage",
        r"\begin{table}[p]",
        r"\centering",
        r"\caption{Bootstrap 95\% coverage (Cov., \%), detection rate (Det., \%),",
        r"and false detection rate (FDR, \%) under the unimodal null (DGP~A6),",
        r"across sample sizes.  Coverage is reported for DGPs with a true",
        r"anti-mode only; ``---'' denotes not applicable.",
        r"Results at $n=250$ are omitted: bootstrap coverage intervals are too",
        r"wide to be informative at that sample size (median half-width",
        r"${>}15$ percentile points).}",
        r"\label{tab:coverage}",
        r"\begin{tabular}{lccccccccc}",
        r"\toprule",
        r" & \multicolumn{3}{c}{$n=500$}"
        r" & \multicolumn{3}{c}{$n=1000$}"
        r" & \multicolumn{3}{c}{$n=5000$}\\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}",
        r"DGP & Cov. & Det. & FDR & Cov. & Det. & FDR & Cov. & Det. & FDR\\",
        r"\midrule",
    ]
    show = {
        "A6": "A6 (unimodal)", "A1": "A1 (weak)",   "A2": "A2 (moderate)",
        "A3": "A3 (sharp)",    "B2": "B2 (heavy tail)",
        "C1": "C1 (joint)",    "C3": "C3 (bathtub)",
    }
    for dgp, label in show.items():
        is_null = (dgp == NULL_DGP)
        row = [label]
        for n in ns_cov:
            sub = df[(df.DGP == dgp) & (df.n == n)]
            if sub.empty:
                row += ["---", "---", "---"]
            elif is_null:
                fdr = _fmt(sub.iloc[0]["FDR"], decimals=1)
                row += ["---", "---", fdr]
            else:
                cov = _fmt(sub.iloc[0]["Coverage"],  decimals=1)
                det = _fmt(sub.iloc[0]["Detection"], decimals=1)
                row += [cov, det, "---"]
        lines.append(" & ".join(row) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TABLES_OUT / "tab_coverage.tex").write_text("\n".join(lines))
    print("  Wrote tab_coverage.tex")


# ── Table 4: Power comparison ─────────────────────────────────────────────────
def write_tab_power(df):
    lines = [
        r"\clearpage",
        r"\begin{table}[p]",
        r"\centering",
        r"\caption{Simulated rejection rates at nominal 5\% level, DGP~A2",
        r"(moderate barrier, $\Delta=2.8$); $R=500$, $B=299$.",
        r"Anti-mode test and dip test applied by the author.",
        r"Z-Dip statistic \citep{ZDip2025} applied by the author;",
        r"Di~Martino et~al.\ validate Z-Dip only on distributions with",
        r"component separation $\geq 8\sigma$ and report no results for",
        r"$\Delta=2.8\sigma$.}",
        r"\label{tab:power}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"$n$ & Anti-mode test & Dip test & Z-Dip\\",
        r"\midrule",
    ]
    sub_a2 = df[df.DGP == "A2"]
    for n in NS:
        row_data = sub_a2[sub_a2.n == n]
        if row_data.empty:
            am = dip = zd = "---"
        else:
            am  = _fmt(row_data.iloc[0]["Power_AM"]   / 100, decimals=3)
            dip = _fmt(row_data.iloc[0]["Power_Dip"]  / 100, decimals=3)
            zd  = _fmt(row_data.iloc[0]["Power_ZDip"] / 100, decimals=3)
        lines.append(f"{n:4d} & {am} & {dip} & {zd}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TABLES_OUT / "tab_power.tex").write_text("\n".join(lines))
    print("  Wrote tab_power.tex")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    missing = [f for f in [
        "table1_bias_rmse.csv", "table2_rates.csv",
        "table3_coverage.csv",  "table4_power.csv",
    ] if not (RESULTS / f).exists()]

    if missing:
        print(f"ERROR: Missing CSV files in results/: {missing}")
        print("Run `python src/simulate.py` first (or `make simulate`).")
        raise SystemExit(1)

    print(f"Reading CSVs from {RESULTS}/ ...")
    df1 = pd.read_csv(RESULTS / "table1_bias_rmse.csv")
    df2 = pd.read_csv(RESULTS / "table2_rates.csv")
    df3 = pd.read_csv(RESULTS / "table3_coverage.csv")
    df4 = pd.read_csv(RESULTS / "table4_power.csv")

    print(f"Writing LaTeX tables to {TABLES_OUT}/ ...")
    write_tab_bias_rmse(df1)
    write_tab_rates(df2)
    write_tab_coverage(df3)
    write_tab_power(df4)
    print("Done. Run `make pdf` to recompile the manuscript.")


if __name__ == "__main__":
    main()
