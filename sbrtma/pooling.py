"""Inverse-variance pooling on the analysis scale (log HR, or a difference).

- dersimonian_laird: fixed effect, DerSimonian-Laird random effects, Q, I2.
- pool_pm_hk: Paule-Mandel tau^2 with a modified Hartung-Knapp interval
  (SE floored at the classic random-effects SE).
- pool_bayes: normal-normal model, half-normal prior on tau, flat prior on
  the mean, marginalised on a grid.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy import stats

Z975 = 1.959963984540054


def t975(df):
    return float(stats.t.ppf(0.975, df))


def _chi2_sf(x, df):
    """Upper-tail chi-square probability (exact)."""
    if df <= 0 or x <= 0:
        return 1.0
    return float(stats.chi2.sf(x, df))


@dataclass
class DLResult:
    k: int
    est_fixed: float
    se_fixed: float
    est_random: float
    se_random: float
    ci_random: tuple
    Q: float
    df: int
    p_Q: float
    I2: float
    tau2: float


def dersimonian_laird(yi, sei):
    yi = np.asarray(yi, float)
    vi = np.asarray(sei, float) ** 2
    k = len(yi)
    if k < 2:
        raise ValueError("need at least 2 studies")
    wi = 1.0 / vi
    est_fixed = np.sum(wi * yi) / np.sum(wi)
    se_fixed = math.sqrt(1.0 / np.sum(wi))
    Q = float(np.sum(wi * (yi - est_fixed) ** 2))
    df = k - 1
    C = np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w = 1.0 / (vi + tau2)
    est_random = np.sum(w * yi) / np.sum(w)
    se_random = math.sqrt(1.0 / np.sum(w))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    ci = (est_random - Z975 * se_random, est_random + Z975 * se_random)
    return DLResult(k, est_fixed, se_fixed, est_random, se_random, ci,
                    Q, df, _chi2_sf(Q, df), I2, tau2)


def leave_one_out(yi, sei):
    """DerSimonian-Laird random-effects estimate omitting each study in turn."""
    yi = np.asarray(yi, float)
    sei = np.asarray(sei, float)
    out = []
    for i in range(len(yi)):
        keep = np.ones(len(yi), bool)
        keep[i] = False
        r = dersimonian_laird(yi[keep], sei[keep])
        out.append((i, r.est_random, r.ci_random, r.I2))
    return out


def paule_mandel(yi, sei, tol=1e-10, maxit=200):
    """tau^2 such that the generalised Q statistic equals k - 1 (bisection)."""
    yi = np.asarray(yi, float)
    vi = np.asarray(sei, float) ** 2
    k = len(yi)
    if k < 2:
        return 0.0
    lo, hi = 0.0, max(10.0, 100 * vi.max())

    def q(t2):
        w = 1.0 / (vi + t2)
        mu = np.sum(w * yi) / np.sum(w)
        return float(np.sum(w * (yi - mu) ** 2))

    if q(0.0) <= k - 1:
        return 0.0
    for _ in range(maxit):
        t2 = 0.5 * (lo + hi)
        if q(t2) > k - 1:
            lo = t2
        else:
            hi = t2
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def pool_pm_hk(yi, sei):
    yi = np.asarray(yi, float)
    sei = np.asarray(sei, float)
    vi, k = sei ** 2, len(yi)
    t2 = paule_mandel(yi, sei)
    w = 1.0 / (vi + t2)
    mu = float(np.sum(w * yi) / np.sum(w))
    se_classic = float(math.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * (yi - mu) ** 2) / (k - 1)) if k > 1 else 1.0
    se_hk = se_classic * math.sqrt(q)
    se = max(se_hk, se_classic)
    t = t975(k - 1)
    t_pi = t975(max(1, k - 2))
    half_pi = t_pi * math.sqrt(t2 + se_classic ** 2)
    return {"k": k, "est": mu, "se": se, "se_classic": se_classic,
            "ci": (mu - t * se, mu + t * se), "tau2": t2,
            "df": k - 1, "t_crit": t, "pred_interval": (mu - half_pi, mu + half_pi),
            "floor_applied": se_hk < se_classic}


def pool_bayes(yi, sei, prior_tau=0.5, n_tau=2000, n_mu=2000):
    yi = np.asarray(yi, float)
    vi = np.asarray(sei, float) ** 2
    taus = np.linspace(1e-6, max(2.0, 8 * prior_tau), n_tau)
    logp = np.empty_like(taus)
    mus = np.empty_like(taus)
    ses = np.empty_like(taus)
    for j, tau in enumerate(taus):
        w = 1.0 / (vi + tau ** 2)
        mu = np.sum(w * yi) / np.sum(w)
        mus[j], ses[j] = mu, math.sqrt(1.0 / np.sum(w))
        # marginal likelihood of tau with mu integrated out under a flat prior
        ll = (-0.5 * np.sum(np.log(vi + tau ** 2))
              - 0.5 * np.sum(w * (yi - mu) ** 2)
              + 0.5 * math.log(1.0 / np.sum(w)))
        logp[j] = ll - 0.5 * (tau / prior_tau) ** 2
    logp -= logp.max()
    p = np.exp(logp)
    p /= np.trapezoid(p, taus)

    grid = np.linspace(min(mus - 6 * ses), max(mus + 6 * ses), n_mu)
    dens = np.zeros_like(grid)
    for j in range(n_tau):
        dens += p[j] * np.exp(-0.5 * ((grid - mus[j]) / ses[j]) ** 2) / (
            ses[j] * math.sqrt(2 * math.pi))
    dens /= np.trapezoid(dens, grid)
    cdf = np.concatenate(([0.0], np.cumsum(
        0.5 * (dens[1:] + dens[:-1]) * (grid[1] - grid[0]))))
    cdf /= cdf[-1]
    qs = np.interp([0.025, 0.5, 0.975], cdf, grid)
    ctau = np.concatenate(([0.0], np.cumsum(
        0.5 * (p[1:] + p[:-1]) * (taus[1] - taus[0]))))
    ctau /= ctau[-1]
    tau_median = float(np.interp(0.5, ctau, taus))
    return {"est": float(qs[1]), "ci": (float(qs[0]), float(qs[2])),
            "tau_median": tau_median, "prior_tau": prior_tau}
