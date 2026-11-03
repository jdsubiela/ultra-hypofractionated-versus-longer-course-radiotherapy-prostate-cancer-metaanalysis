"""Reconstruction of individual patient data from a digitized Kaplan-Meier
curve and its numbers-at-risk table.

Python port of the iterative algorithm of Guyot P, Ades AE, Ouwens MJ,
Welton NJ. BMC Med Res Methodol 2012;12:9.
"""
import numpy as np


def prepare_curve(t, s, maxy):
    """Sort, rescale to 0-1, enforce a non-increasing curve starting at (0, 1)."""
    t = np.asarray(t, float)
    s = np.asarray(s, float)
    if maxy != 1.0:
        s = s / maxy
    order = np.argsort(t)
    t, s = t[order], np.minimum.accumulate(s[order])
    s = np.clip(s, 0, 1)
    if t[0] > 0:
        t = np.insert(t, 0, 0.0)
        s = np.insert(s, 0, 1.0)
    s[0] = 1.0
    return t, s


def _spread_censor(t_S, lo_idx, hi_idx, m):
    """Counts per bin of m censoring times spaced evenly between two points."""
    if m <= 0:
        return np.zeros(hi_idx - lo_idx)
    cen_t = t_S[lo_idx] + np.arange(1, m + 1) * (t_S[hi_idx] - t_S[lo_idx]) / (m + 1)
    counts, _ = np.histogram(cen_t, bins=t_S[lo_idx:hi_idx + 1])
    return counts


def reconstruct_ipd(t_S, S, t_risk, n_risk, tot_events=None):
    """Return an (N, 2) array of [time, status] (status 1 = event, 0 = censored).

    t_S, S: digitized curve, sorted by time, S in 0-1 and non-increasing.
    t_risk, n_risk: numbers-at-risk table.
    tot_events: total events in the arm, or None.
    """
    t_S = np.asarray(t_S, float)
    S = np.asarray(S, float)
    t_risk = np.asarray(t_risk, float)
    n_risk = np.array(n_risk, float)
    n = len(t_S)
    n_int = len(n_risk)

    # first and last digitized index inside each at-risk interval
    lower = np.zeros(n_int, dtype=int)
    upper = np.zeros(n_int, dtype=int)
    for i in range(n_int):
        lower[i] = int(np.searchsorted(t_S, t_risk[i], side="left"))
        if i < n_int - 1:
            upper[i] = int(np.searchsorted(t_S, t_risk[i + 1], side="left")) - 1
        else:
            upper[i] = n - 1
    lower = np.clip(lower, 0, n - 1)
    upper = np.clip(upper, 0, n - 1)

    n_censor = np.zeros(n_int)
    n_hat = np.full(n + 1, n_risk[0] + 1)
    cen = np.zeros(n)
    d = np.zeros(n)
    KM_hat = np.ones(n)
    last_i = np.full(n_int, lower[0], dtype=int)

    def step(lo, hi, last, first_is_origin, clamp):
        for k in range(lo, hi + 1):
            if first_is_origin and k == lo:
                d[k] = 0
                KM_hat[k] = 1.0
            else:
                d[k] = round(n_hat[k] * (1 - (S[k] / KM_hat[last])))
            KM_hat[k] = KM_hat[last] * (1 - (d[k] / n_hat[k])) if n_hat[k] > 0 else KM_hat[last]
            n_hat[k + 1] = n_hat[k] - d[k] - cen[k]
            if clamp and n_hat[k + 1] < 0:
                n_hat[k + 1] = 0
                cen[k] = n_hat[k] - d[k]
            if d[k] != 0:
                last = k
        return last

    # all intervals except the last
    for i in range(n_int - 1):
        if KM_hat[lower[i]] != 0:
            n_censor[i] = round(n_risk[i] * S[lower[i + 1]] / KM_hat[lower[i]] - n_risk[i + 1])
        else:
            n_censor[i] = 0
        for _ in range(10000):
            if n_censor[i] <= 0:
                cen[lower[i]:upper[i] + 1] = 0
                n_censor[i] = 0
            else:
                counts = _spread_censor(t_S, lower[i], lower[i + 1], int(n_censor[i]))
                cen[lower[i]:upper[i] + 1] = 0
                cen[lower[i]:lower[i + 1]] = counts
            n_hat[lower[i]] = n_risk[i]
            last = step(lower[i], upper[i], last_i[i], i == 0, clamp=False)
            nh = n_hat[lower[i + 1]]
            if nh > n_risk[i + 1] or (nh < n_risk[i + 1] and n_censor[i] > 0):
                n_censor[i] = n_censor[i] + (nh - n_risk[i + 1])
            else:
                break
        if n_hat[lower[i + 1]] < n_risk[i + 1]:
            n_risk[i + 1] = n_hat[lower[i + 1]]
        last_i[i + 1] = last

    # last interval: censoring rate carried over from the earlier intervals
    ni = n_int - 1
    if n_int > 1:
        span = t_S[lower[ni]] - t_S[lower[0]]
        rate = (np.sum(n_censor[:ni]) * (t_S[upper[ni]] - t_S[lower[ni]]) / span) if span > 0 else 0
        n_censor[ni] = min(round(rate), n_risk[ni])
    else:
        n_censor[ni] = 0

    def fill_last_censor():
        if n_censor[ni] <= 0:
            cen[lower[ni]:upper[ni]] = 0
            n_censor[ni] = 0
        else:
            cen[lower[ni]:upper[ni]] = _spread_censor(t_S, lower[ni], upper[ni],
                                                      int(n_censor[ni]))

    fill_last_censor()
    n_hat[lower[ni]] = n_risk[ni]
    step(lower[ni], upper[ni], last_i[ni], n_int == 1, clamp=True)

    # match the reported total number of events
    if tot_events is not None:
        sumdL = 0
        if n_int > 1:
            sumdL = np.sum(d[:upper[ni - 1] + 1])
            if sumdL >= tot_events:
                d[lower[ni]:upper[ni] + 1] = 0
                cen[lower[ni]:upper[ni]] = 0
                n_hat[lower[ni] + 1:upper[ni] + 2] = n_risk[ni]
        if sumdL < tot_events or n_int == 1:
            sumd = np.sum(d[:upper[ni] + 1])
            for _ in range(10000):
                if not (sumd > tot_events or (sumd < tot_events and n_censor[ni] > 0)):
                    break
                n_censor[ni] = n_censor[ni] + (sumd - tot_events)
                fill_last_censor()
                n_hat[lower[ni]] = n_risk[ni]
                step(lower[ni], upper[ni], last_i[ni], n_int == 1, clamp=True)
                sumd = np.sum(d[:upper[ni] + 1])

    rows = []
    for k in range(n):
        rows += [(t_S[k], 1)] * int(d[k])
        rows += [(t_S[k], 0)] * int(cen[k])
    # patients still at risk after the last point are censored there
    remaining = int(round(n_hat[upper[n_int - 1] + 1]))
    if remaining > 0:
        rows += [(t_S[upper[n_int - 1]], 0)] * remaining
    return np.array(rows, dtype=float) if rows else np.empty((0, 2))


def km_step(time, status):
    """Kaplan-Meier step function: unique times and survival."""
    order = np.argsort(time)
    t, s = np.asarray(time)[order], np.asarray(status)[order]
    uniq = np.unique(t)
    surv, S = [], 1.0
    for ut in uniq:
        at_risk = (t >= ut).sum()
        d = ((t == ut) & (s == 1)).sum()
        if at_risk > 0 and d > 0:
            S *= (1 - d / at_risk)
        surv.append(S)
    return uniq, np.array(surv)


def km_at(uniq, surv, q):
    idx = np.searchsorted(uniq, q, side="right") - 1
    return 1.0 if idx < 0 else surv[idx]
