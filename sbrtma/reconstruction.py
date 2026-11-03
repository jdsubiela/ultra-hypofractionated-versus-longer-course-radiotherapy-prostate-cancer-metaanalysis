"""Reconstruct IPD for every arm and run the validation checks.

Checks per arm:
  - maximum distance between the reconstructed KM and the digitized points
    (flag if > 0.03);
  - reconstructed vs reported total events, when the total is used
    (flag if off by more than max(3, 5%));
  - digitized and reconstructed survival at each printed milestone vs the
    published value (flag if either is off by more than 0.5 percentage points).
A pair with a flagged arm is excluded from the IPD analyses.
"""
import csv
import os

import numpy as np

from .guyot import km_at, km_step, prepare_curve, reconstruct_ipd

FIT_TOL = 0.03
MILESTONE_TOL = 0.5


def read_points(path):
    with open(path) as fh:
        return [tuple(map(float, ln.split(","))) for ln in fh if ln.strip()]


def reconstruct_arm(cid, c, digitized_dir, out_dir):
    pts = read_points(os.path.join(digitized_dir, cid + ".csv"))
    t_risk, n_risk = c["t_risk"], c["n_risk"]
    tmax = c.get("t_max")
    if tmax is not None:
        keep = [i for i, tv in enumerate(t_risk) if tv <= tmax + 1e-9]
        t_risk = [t_risk[i] for i in keep]
        n_risk = [n_risk[i] for i in keep]
        pts = [p for p in pts if p[0] <= tmax + 1e-9]
    tot_events = c.get("events") if tmax is None else None

    t_S, S = prepare_curve([p[0] for p in pts], [p[1] for p in pts], c["maxy"])
    ipd = reconstruct_ipd(t_S, S, t_risk, n_risk, tot_events=tot_events)

    with open(os.path.join(out_dir, cid + "_ipd.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "status"])
        w.writerows(ipd.tolist())

    uniq, surv = km_step(ipd[:, 0], ipd[:, 1])
    rec_events = int(ipd[:, 1].sum())
    fit_dev = float(max(abs(km_at(uniq, surv, t) - s) for t, s in zip(t_S, S)))
    flags = []
    if fit_dev > FIT_TOL:
        flags.append("curve fit")
    if tot_events and abs(rec_events - tot_events) > max(3, 0.05 * tot_events):
        flags.append("event count")

    milestones = []
    for key, yr in (("reported_3y", 3), ("reported_5y", 5), ("reported_10y", 10)):
        if c.get(key) is None or (tmax is not None and yr > tmax + 1e-9):
            continue
        v = [s for t, s in pts if t <= yr + 1e-9]
        dig = v[-1] if v else 100.0
        rec = float(km_at(uniq, surv, yr)) * 100.0
        ok = abs(dig - c[key]) <= MILESTONE_TOL and abs(rec - c[key]) <= MILESTONE_TOL
        milestones.append({"year": yr, "published": c[key],
                           "digitized": round(dig, 2), "diff": round(dig - c[key], 2),
                           "reconstructed": round(rec, 2),
                           "diff_reconstructed": round(rec - c[key], 2), "pass": ok})
        if not ok:
            flags.append(f"milestone {yr} y")

    return {"curve_id": cid, "trial": c["trial"], "arm": c["arm"],
            "outcome": c["outcome"], "rows": int(len(ipd)),
            "n_at_risk_t0": int(n_risk[0]), "events_reconstructed": rec_events,
            "events_reported": tot_events, "max_fit_deviation": round(fit_dev, 4),
            "milestones": milestones, "flags": flags, "pass": not flags}


def reconstruct_all(curves, digitized_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    return {cid: reconstruct_arm(cid, c, digitized_dir, out_dir)
            for cid, c in curves.items()}
