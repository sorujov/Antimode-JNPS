#!/usr/bin/env python3
"""
Empirical illustration: Old Faithful geyser eruption durations.
Applies the anti-mode detection, bootstrap CI, and both tests
(sup-norm and pointwise) to the classic bimodal dataset.

Data: 272 eruption durations (minutes) from Azzalini & Bowman (1990).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings("ignore")

# ── Old Faithful eruption durations (minutes) ──
# Source: R datasets::faithful$eruptions (272 observations)
# Azzalini, A. and Bowman, A. W. (1990). A look at some data on the
# Old Faithful geyser. Applied Statistics, 39, 357–365.
FAITHFUL = np.array([
    3.600, 1.800, 3.333, 2.283, 4.533, 2.883, 4.700, 3.600, 1.950,
    4.350, 1.833, 3.917, 4.200, 1.750, 4.700, 2.167, 1.750, 4.800,
    1.600, 4.250, 1.800, 1.750, 3.450, 3.067, 4.533, 3.600, 1.967,
    4.083, 3.850, 4.433, 4.300, 4.467, 3.367, 4.033, 3.833, 2.017,
    1.867, 4.833, 1.833, 4.783, 4.350, 1.883, 4.567, 1.750, 4.533,
    3.317, 3.833, 2.100, 4.633, 2.000, 4.800, 4.716, 1.833, 4.833,
    1.733, 4.883, 3.717, 1.667, 4.567, 4.317, 2.233, 4.500, 1.750,
    4.800, 1.817, 4.400, 4.167, 4.700, 2.067, 4.700, 4.033, 1.967,
    4.500, 4.000, 1.983, 5.067, 2.017, 4.567, 3.883, 3.600, 4.133,
    4.333, 4.100, 2.633, 4.067, 4.933, 3.950, 4.517, 2.167, 4.000,
    2.200, 4.333, 1.867, 4.817, 1.833, 4.300, 4.667, 3.750, 1.867,
    4.900, 2.483, 4.367, 2.100, 4.500, 4.050, 1.867, 4.700, 1.783,
    4.850, 3.683, 4.733, 2.300, 4.900, 4.417, 1.700, 4.633, 2.317,
    4.600, 1.817, 4.417, 2.617, 4.067, 4.250, 1.967, 4.600, 3.767,
    1.917, 4.500, 2.267, 4.650, 1.867, 4.167, 2.800, 4.333, 1.833,
    4.383, 1.883, 4.933, 2.033, 3.733, 4.233, 2.233, 4.533, 4.817,
    4.333, 1.983, 4.633, 2.017, 5.100, 1.800, 5.033, 4.000, 2.400,
    4.600, 3.567, 4.000, 4.500, 4.083, 1.800, 3.967, 2.200, 4.150,
    2.000, 3.833, 3.500, 4.583, 2.367, 5.000, 1.933, 4.617, 1.917,
    2.083, 4.583, 3.833, 4.167, 4.333, 4.500, 2.417, 4.000, 4.167,
    1.883, 4.583, 4.250, 3.767, 2.033, 4.433, 4.083, 1.833, 4.417,
    2.183, 4.800, 1.833, 4.800, 4.100, 3.966, 4.233, 3.500, 4.366,
    2.250, 4.667, 2.100, 4.350, 4.133, 1.867, 4.600, 1.783, 4.367,
    3.850, 1.933, 4.500, 2.383, 4.700, 1.867, 3.833, 3.417, 4.233,
    2.400, 4.800, 2.000, 4.150, 1.867, 4.267, 1.750, 4.483, 4.000,
    4.117, 4.083, 4.267, 3.917, 4.550, 4.083, 2.417, 4.183, 2.217,
    4.450, 1.883, 1.850, 4.283, 3.950, 2.333, 4.150, 2.350, 4.933,
    2.900, 4.583, 3.833, 2.083, 4.367, 2.133, 4.350, 2.200, 4.450,
    3.567, 4.500, 4.150, 3.817, 3.917, 4.450, 2.000, 4.283, 4.767,
    4.533, 1.850, 4.250, 1.983, 2.250, 4.750, 4.117, 2.150, 4.417,
    1.817, 4.467,
])

# ── Parameters (same as main simulation) ──
DELTA = 0.05
TAU   = 0.05
KAPPA = 0.05
B     = 999   # more bootstrap reps for final results
SEED  = 42


def detect(x):
    n  = len(x)
    h  = 0.6 * 1.06 * x.std() * n**(-0.2)
    pg = np.linspace(DELTA, 1 - DELTA, 600)
    qg = np.quantile(x, pg)
    f  = np.clip(gaussian_kde(x, bw_method=h / x.std()).evaluate(qg), 1e-10, None)
    qd = 1.0 / f
    pk = [i for i in range(1, len(qd) - 1)
          if qd[i] > qd[i-1] and qd[i] > qd[i+1]]
    if not pk:
        return [], pg, qd
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
    return [p for p, _ in merged], pg, qd


def bci(x, ph, rng):
    bt = []
    for _ in range(B):
        xb = rng.choice(x, len(x), replace=True)
        a, _, _ = detect(xb)
        if a:
            bt.append(min(a, key=lambda p: abs(p - ph)))
    if len(bt) < 5:
        return np.nan, np.nan
    return np.percentile(bt, 2.5), np.percentile(bt, 97.5)


def am_test_sup(x, rng):
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
        xb     = rng.normal(mu0, s0, n)
        mb, sb = xb.mean(), xb.std()
        hb     = 0.6 * 1.06 * sb * n**(-0.2)
        fb     = np.clip(gaussian_kde(xb, bw_method=hb / sb).evaluate(np.quantile(xb, pg)), 1e-10, None)
        qdb    = 1.0 / fb
        qnb    = sb / stats.norm.pdf(stats.norm.ppf(pg))
        null.append(np.max(np.clip(qdb - qnb, 0, None)))
    return T > np.percentile(null, 95)


def am_test_pointwise(x, rng):
    n   = len(x)
    pg  = np.linspace(DELTA, 1 - DELTA, 300)
    h   = 0.6 * 1.06 * x.std() * n**(-0.2)
    qg  = np.quantile(x, pg)
    f   = np.clip(gaussian_kde(x, bw_method=h / x.std()).evaluate(qg), 1e-10, None)
    qd  = 1.0 / f
    mu0, s0 = x.mean(), x.std()
    qn  = s0 / stats.norm.pdf(stats.norm.ppf(pg))
    diff = qd - qn
    idx_star = np.argmax(diff)
    T = max(diff[idx_star], 0.0)
    if T == 0.0:
        return False
    null = []
    for _ in range(B):
        xb     = rng.normal(mu0, s0, n)
        mb, sb = xb.mean(), xb.std()
        hb     = 0.6 * 1.06 * sb * n**(-0.2)
        qgb    = np.quantile(xb, pg)
        fb     = np.clip(gaussian_kde(xb, bw_method=hb / sb).evaluate(qgb), 1e-10, None)
        qdb    = 1.0 / fb
        qnb    = sb / stats.norm.pdf(stats.norm.ppf(pg))
        null.append(max(qdb[idx_star] - qnb[idx_star], 0.0))
    return T > np.percentile(null, 95)


def dip_test(x):
    import diptest
    _, pval = diptest.diptest(x)
    return pval


if __name__ == "__main__":
    x = FAITHFUL
    n = len(x)
    rng = np.random.default_rng(SEED)

    print("=" * 60)
    print("EMPIRICAL ILLUSTRATION: Old Faithful Eruption Durations")
    print("=" * 60)
    print(f"n = {n}")
    print(f"Mean = {x.mean():.3f}, Std = {x.std():.3f}")
    print(f"Min = {x.min():.3f}, Max = {x.max():.3f}")
    print()

    # Detection
    am, pg, qd = detect(x)
    print(f"Detected anti-modes: {len(am)}")
    for i, p in enumerate(am):
        xval = np.quantile(x, p)
        print(f"  Anti-mode {i+1}: p* = {p:.4f}, x* = Q(p*) = {xval:.3f} minutes")

    # Bootstrap CI
    if am:
        ph = am[0]
        lo, hi = bci(x, ph, np.random.default_rng(SEED + 1))
        xlo = np.quantile(x, lo) if not np.isnan(lo) else np.nan
        xhi = np.quantile(x, hi) if not np.isnan(hi) else np.nan
        print(f"\n  95% Bootstrap CI for p*: [{lo:.4f}, {hi:.4f}]")
        print(f"  95% Bootstrap CI for x*: [{xlo:.3f}, {xhi:.3f}] minutes")

    # Tests
    print("\nHypothesis tests (H0: unimodal):")
    rej_sup = am_test_sup(x, np.random.default_rng(SEED + 2))
    rej_pw  = am_test_pointwise(x, np.random.default_rng(SEED + 3))
    dip_p   = dip_test(x)
    print(f"  Sup-norm test:    {'REJECT' if rej_sup else 'fail to reject'} at 5%")
    print(f"  Pointwise test:   {'REJECT' if rej_pw  else 'fail to reject'} at 5%")
    print(f"  Dip test:         p-value = {dip_p:.4f} ({'REJECT' if dip_p < 0.05 else 'fail to reject'} at 5%)")

    # Fisher information profile at anti-mode
    if am:
        h = 0.6 * 1.06 * x.std() * n**(-0.2)
        xstar = np.quantile(x, am[0])
        fstar = gaussian_kde(x, bw_method=h / x.std()).evaluate([xstar])[0]
        Istar = fstar**2 / (am[0] * (1 - am[0]))
        print(f"\n  Fisher info profile at anti-mode: I(p*) = {Istar:.4f}")
        print(f"  Density at anti-mode: f(x*) = {fstar:.4f}")
        print(f"  Quantile density at anti-mode: q(p*) = {1/fstar:.4f}")

    print("\nDone.")
