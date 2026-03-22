#!/usr/bin/env python3
"""Quick prototype: pointwise anti-mode test vs sup-norm test on DGP A2.
Evaluates T_n only at the detected anti-mode percentile p_hat* instead of
taking the supremum over all p. Should have much higher power at moderate n.
"""
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde

DELTA = 0.05
TAU   = 0.05
KAPPA = 0.05
B     = 299
R     = 50
NS    = [250, 500, 1000, 5000]
SEED  = 42


def sample_A2(n, rng):
    z = rng.binomial(1, .5, n)
    return np.where(z, rng.normal(0, 1, n), rng.normal(2.8, 1, n))


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


def am_test_sup(x, rng):
    """Original sup-norm test."""
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
    """Pointwise test: bootstrap the full procedure under H0.
    For each bootstrap sample from the fitted null, detect the argmax
    of (qd_boot - qn_boot) and record that value. This accounts for
    the selection bias of choosing which percentile to evaluate."""
    n   = len(x)
    pg  = np.linspace(DELTA, 1 - DELTA, 300)
    h   = 0.6 * 1.06 * x.std() * n**(-0.2)
    qg  = np.quantile(x, pg)
    f   = np.clip(gaussian_kde(x, bw_method=h / x.std()).evaluate(qg), 1e-10, None)
    qd  = 1.0 / f
    mu0, s0 = x.mean(), x.std()
    qn  = s0 / stats.norm.pdf(stats.norm.ppf(pg))

    # Observed: max pointwise difference
    diff = qd - qn
    idx_star = np.argmax(diff)
    T = max(diff[idx_star], 0.0)

    if T == 0.0:
        return False

    # Bootstrap null: for each null sample, find ITS OWN argmax of (qdb - qnb)
    null = []
    for _ in range(B):
        xb     = rng.normal(mu0, s0, n)
        mb, sb = xb.mean(), xb.std()
        hb     = 0.6 * 1.06 * sb * n**(-0.2)
        qgb    = np.quantile(xb, pg)
        fb     = np.clip(gaussian_kde(xb, bw_method=hb / sb).evaluate(qgb), 1e-10, None)
        qdb    = 1.0 / fb
        qnb    = sb / stats.norm.pdf(stats.norm.ppf(pg))
        diff_b = qdb - qnb
        null.append(max(diff_b[np.argmax(diff_b)], 0.0))
    return T > np.percentile(null, 95)


def _one_rep_A2(args):
    n, r = args
    rng   = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13)
    x     = sample_A2(n, rng)
    rng_s = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13 + 1)
    rng_p = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13 + 2)
    return int(am_test_sup(x, rng_s)), int(am_test_pointwise(x, rng_p))


def _one_rep_null(args):
    """Size check under unimodal null N(0,1)."""
    n, r = args
    rng   = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13 + 100)
    x     = rng.normal(0, 1, n)
    rng_s = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13 + 101)
    rng_p = np.random.default_rng(SEED + r * 997 + NS.index(n) * 13 + 102)
    return int(am_test_sup(x, rng_s)), int(am_test_pointwise(x, rng_p))


if __name__ == "__main__":
    import os
    import warnings
    from multiprocessing import Pool
    warnings.filterwarnings("ignore")

    NCPUS = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))

    print(f"Prototype: pointwise vs sup-norm test")
    print(f"R={R}, B={B}, workers={NCPUS}")

    with Pool(processes=NCPUS) as pool:
        # Power under A2
        print(f"\n=== POWER (DGP A2, Delta=2.8) ===")
        print(f"{'n':>6}  {'Sup-norm':>10}  {'Pointwise':>10}")
        print("-" * 32)
        for n in NS:
            args = [(n, r) for r in range(R)]
            results = pool.map(_one_rep_A2, args)
            rej_sup = sum(r[0] for r in results)
            rej_pw  = sum(r[1] for r in results)
            print(f"{n:>6}  {rej_sup/R:>10.3f}  {rej_pw/R:>10.3f}")

        # Size under null N(0,1)
        print(f"\n=== SIZE (DGP A6, unimodal N(0,1)) ===")
        print(f"{'n':>6}  {'Sup-norm':>10}  {'Pointwise':>10}")
        print("-" * 32)
        for n in NS:
            args = [(n, r) for r in range(R)]
            results = pool.map(_one_rep_null, args)
            rej_sup = sum(r[0] for r in results)
            rej_pw  = sum(r[1] for r in results)
            print(f"{n:>6}  {rej_sup/R:>10.3f}  {rej_pw/R:>10.3f}")

    print("\nDone.")
