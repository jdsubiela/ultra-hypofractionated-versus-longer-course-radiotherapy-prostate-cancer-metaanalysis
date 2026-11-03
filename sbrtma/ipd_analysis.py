"""Analyses on the reconstructed IPD.

Per pair: Cox HR over full follow-up and with follow-up truncated at 3 and
5 years (Grambsch-Therneau test of proportional hazards), and RMST at 3, 5
and 10 years (Klein-Moeschberger variance). Across trials: common-effect
pooling of the truncated HRs (with Paule-Mandel/Hartung-Knapp alongside) and
DerSimonian-Laird random-effects pooling of the RMST differences.
"""
import math
import os

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

from . import pooling

LANDMARKS = [3.0, 5.0]
RMST_TAUS = [3.0, 5.0, 10.0]
# A horizon is analysed if reconstructed follow-up ends within this many years
# of it (digitized time axes do not land exactly on the horizon).
LANDMARK_TOL = 0.02

GROUPS = {
    "composite_failure": ["Failure-free survival", "Disease-free survival",
                          "Freedom from biochemical/clinical failure"],
    "overall_survival": ["Overall survival"],
}


def load_ipd(ipd_dir, cid):
    df = pd.read_csv(os.path.join(ipd_dir, cid + "_ipd.csv"))
    return df["time"].to_numpy(float), df["status"].to_numpy(float)


def km_table(t, d):
    order = np.argsort(t)
    t, d = t[order], d[order]
    out_t, out_s, out_n, out_d = [], [], [], []
    S = 1.0
    for ut in np.unique(t):
        n_i = int((t >= ut).sum())
        d_i = int(((t == ut) & (d == 1)).sum())
        if d_i:
            S *= (1 - d_i / n_i)
        out_t.append(float(ut)); out_s.append(S); out_n.append(n_i); out_d.append(d_i)
    return np.array(out_t), np.array(out_s), np.array(out_n), np.array(out_d)


def rmst(t, d, tau):
    """RMST up to tau and its SE (Klein-Moeschberger)."""
    ts, S, n, ev = km_table(t, d)
    keep = ts <= tau
    ts, S, n, ev = ts[keep], S[keep], n[keep], ev[keep]
    area = float(np.sum(np.diff(np.concatenate(([0.0], ts, [tau])))
                        * np.concatenate(([1.0], S))))
    var = 0.0
    for i in range(len(ts)):
        if ev[i] == 0 or n[i] <= ev[i]:
            continue
        tail = float(np.sum(np.diff(np.concatenate((ts[i:], [tau]))) * S[i:]))
        var += tail ** 2 * ev[i] / (n[i] * (n[i] - ev[i]))
    return area, math.sqrt(max(var, 0.0))


def cox(t_e, d_e, t_c, d_c, tau=None):
    """HR experimental vs control; tau truncates follow-up at that time."""
    te, de, tc, dc = t_e.copy(), d_e.copy(), t_c.copy(), d_c.copy()
    if tau is not None:
        de = np.where(te > tau, 0.0, de); te = np.minimum(te, tau)
        dc = np.where(tc > tau, 0.0, dc); tc = np.minimum(tc, tau)
    df = pd.DataFrame({"T": np.concatenate([te, tc]),
                       "E": np.concatenate([de, dc]),
                       "arm": np.concatenate([np.ones(len(te)), np.zeros(len(tc))])})
    cph = CoxPHFitter()
    cph.fit(df, duration_col="T", event_col="E")
    hr = float(np.exp(cph.params_["arm"]))
    se = float(cph.standard_errors_["arm"])
    lo, hi = math.exp(math.log(hr) - 1.959964 * se), math.exp(math.log(hr) + 1.959964 * se)
    p = float(cph.summary.loc["arm", "p"])
    ph = proportional_hazard_test(cph, df, time_transform="km")
    p_ph = float(np.asarray(ph.p_value).ravel()[0])
    return {"hr": round(hr, 4), "se_log": round(se, 5),
            "lo": round(lo, 4), "hi": round(hi, 4), "p": round(p, 5),
            "p_proportional_hazards": round(p_ph, 5),
            "n": int(len(df)), "events": int(df["E"].sum())}


def analyse_pair(te, de, tc, dc):
    tmax = min(te.max(), tc.max())
    entry = {"t_max": round(float(tmax), 2),
             "n_exp": int(len(te)), "n_ctl": int(len(tc)),
             "events_exp": int(de.sum()), "events_ctl": int(dc.sum()),
             "rmst": [], "cox": {"full": cox(te, de, tc, dc)}}
    for lm in LANDMARKS:
        if lm > tmax + LANDMARK_TOL:
            continue
        c = cox(te, de, tc, dc, tau=lm)
        if lm > tmax:
            c["landmark_tolerance_y"] = round(float(lm - tmax), 5)
        entry["cox"][f"landmark_{lm:g}y"] = c
    for tau in RMST_TAUS:
        if tau > tmax + LANDMARK_TOL:
            continue
        a_e, s_e = rmst(te, de, tau)
        a_c, s_c = rmst(tc, dc, tau)
        diff = a_e - a_c
        se = math.sqrt(s_e ** 2 + s_c ** 2)
        lo, hi = diff - 1.959964 * se, diff + 1.959964 * se
        z = diff / se if se > 0 else 0.0
        p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
        entry["rmst"].append({
            "tau_y": tau, "rmst_exp": round(a_e, 4), "rmst_ctl": round(a_c, 4),
            "diff_y": round(diff, 4), "se": round(se, 4),
            "lo": round(lo, 4), "hi": round(hi, 4), "p": round(p, 5),
            "diff_months": round(diff * 12, 2),
            "lo_months": round(lo * 12, 2), "hi_months": round(hi * 12, 2)})
    return entry


def pool_landmark(pairs, outcomes, lm):
    yi, sei, trials, tolerated = [], [], [], []
    for e in pairs:
        if e["excluded"] or e["outcome"] not in outcomes:
            continue
        c = e["cox"].get(f"landmark_{lm:g}y")
        if not c:
            continue
        yi.append(math.log(c["hr"])); sei.append(c["se_log"]); trials.append(e["trial"])
        if "landmark_tolerance_y" in c:
            tolerated.append(e["trial"])
    if len(yi) < 2:
        return None
    yi, sei = np.array(yi), np.array(sei)
    dl = pooling.dersimonian_laird(yi, sei)
    w = 1.0 / sei ** 2
    mu = float(np.sum(w * yi) / np.sum(w))
    se = math.sqrt(1.0 / np.sum(w))
    pm = pooling.pool_pm_hk(yi, sei)
    out = {"k": dl.k, "trials": trials, "model": "common effect",
           "hr": round(math.exp(mu), 3),
           "lo": round(math.exp(mu - 1.96 * se), 3),
           "hi": round(math.exp(mu + 1.96 * se), 3),
           "I2": round(dl.I2, 1), "tau2": round(dl.tau2, 5), "p_Q": round(dl.p_Q, 4),
           "random_pm_hk": {"hr": round(math.exp(pm["est"]), 3),
                            "lo": round(math.exp(pm["ci"][0]), 3),
                            "hi": round(math.exp(pm["ci"][1]), 3),
                            "tau2": round(pm["tau2"], 5)}}
    if tolerated:
        out["landmark_tolerance_trials"] = tolerated
    return out


def pool_rmst(pairs, outcomes, tau):
    yi, sei, trials = [], [], []
    for e in pairs:
        if e["excluded"] or e["outcome"] not in outcomes:
            continue
        row = next((r for r in e["rmst"] if r["tau_y"] == tau), None)
        if row:
            yi.append(row["diff_y"]); sei.append(row["se"]); trials.append(e["trial"])
    if len(yi) < 2:
        return None
    r = pooling.dersimonian_laird(np.array(yi), np.array(sei))
    return {"k": r.k, "trials": trials,
            "diff_months": round(r.est_random * 12, 2),
            "lo_months": round(r.ci_random[0] * 12, 2),
            "hi_months": round(r.ci_random[1] * 12, 2),
            "I2": round(r.I2, 1), "tau2": round(r.tau2, 6)}


def run(pairs_def, qa, ipd_dir):
    pairs = []
    for trial, outcome, exp_id, ctl_id in pairs_def:
        failed = [c for c in (exp_id, ctl_id) if not qa[c]["pass"]]
        if failed:
            pairs.append({"trial": trial, "outcome": outcome, "excluded": True,
                          "failed_validation": failed})
            continue
        te, de = load_ipd(ipd_dir, exp_id)
        tc, dc = load_ipd(ipd_dir, ctl_id)
        pairs.append({"trial": trial, "outcome": outcome, "excluded": False,
                      **analyse_pair(te, de, tc, dc)})

    pooled = {}
    for gname, outcomes in GROUPS.items():
        pooled[gname] = {}
        for lm in LANDMARKS:
            r = pool_landmark(pairs, outcomes, lm)
            if r:
                pooled[gname][f"landmark_{lm:g}y"] = r
        for tau in RMST_TAUS:
            r = pool_rmst(pairs, outcomes, tau)
            if r:
                pooled[gname][f"rmst_{tau:g}y"] = r
    return {"pairs": pairs, "meta": pooled}
