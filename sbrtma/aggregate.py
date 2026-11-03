"""Aggregate-data synthesis from the extracted trial results.

Primary model: Paule-Mandel tau^2 with a modified Hartung-Knapp interval.
Sensitivity: DerSimonian-Laird and a Bayesian normal-normal model.
"""
import math

import numpy as np
from scipy import stats

from . import pooling

POOLED_OUTCOMES = ("efficacy_primary", "biochemical_failure", "late_gu_g2",
                   "late_gi_g2", "overall_survival")


def eqd2(total_dose, n_fx, ab):
    """EQD2 = D * (d + a/b) / (2 + a/b), with d = D / n."""
    d = total_dose / n_fx
    return total_dose * (d + ab) / (2.0 + ab)


def _eqd2_row(study, arm, dose, nfx):
    return {"study": study, "arm": arm, "total_dose": dose, "n_fx": nfx,
            "dose_per_fx": round(dose / nfx, 3),
            "eqd2_ab1_5": round(eqd2(dose, nfx, 1.5), 1),
            "eqd2_ab3": round(eqd2(dose, nfx, 3.0), 1)}


def eqd2_table(studies):
    rows = []
    for key, s in studies.items():
        rows.append(_eqd2_row(key, "experimental", s["dose_exp_gy"], s["nfx_exp"]))
        rows.append(_eqd2_row(key, "control", s["dose_ctrl_gy"], s["nfx_ctrl"]))
        if s.get("internal_boost"):
            rows.append(_eqd2_row(key, "experimental (CTV boost)",
                                  s["internal_boost"], s["nfx_exp"]))
        if s.get("dose_cap_gy"):
            rows.append(_eqd2_row(key, "experimental (PTV max cap)",
                                  s["dose_cap_gy"], s["nfx_exp"]))
    return rows


def log_hr(d):
    """ln(HR) and its SE from the reported 95% CI."""
    if d["ci"] != 95:
        raise ValueError("only 95% confidence intervals are used")
    return math.log(d["hr"]), (math.log(d["hi"]) - math.log(d["lo"])) / (2 * pooling.Z975)


def _exp_ci(ci):
    return [round(math.exp(ci[0]), 3), round(math.exp(ci[1]), 3)]


def pool_outcome(key, outcome, studies):
    names = list(outcome["data"])
    ys = [log_hr(outcome["data"][n]) for n in names]
    yi = np.array([y for y, _ in ys])
    sei = np.array([s for _, s in ys])
    dl = pooling.dersimonian_laird(yi, sei)
    pm = pooling.pool_pm_hk(yi, sei)
    bay = pooling.pool_bayes(yi, sei)
    w_pm = 1.0 / (sei ** 2 + pm["tau2"])

    per_study = []
    for i, name in enumerate(names):
        d = outcome["data"][name]
        per_study.append({
            "study": name, "label": studies[name]["label"],
            "hr": d["hr"], "lo": d["lo"], "hi": d["hi"],
            "yi": float(yi[i]), "sei": float(sei[i]),
            "weight_pct": round(float(w_pm[i] / w_pm.sum() * 100), 1),
            "comparator_type": studies[name]["comparator_type"],
            "src": d["src"]})

    loo = []
    if len(names) >= 3:
        for i, est, ci, i2 in pooling.leave_one_out(yi, sei):
            loo.append({"omitted": names[i], "hr": round(math.exp(est), 3),
                        "lo": round(math.exp(ci[0]), 3), "hi": round(math.exp(ci[1]), 3),
                        "I2": round(i2, 1)})

    return {
        "outcome": key, "label": outcome["label"], "k": dl.k,
        "estimator": "Paule-Mandel + modified Hartung-Knapp",
        "hr_random": round(math.exp(pm["est"]), 3),
        "ci_random": _exp_ci(pm["ci"]),
        "hr_fixed": round(math.exp(dl.est_fixed), 3),
        "pred_interval": _exp_ci(pm["pred_interval"]),
        "pred_interval_informative": dl.k >= 5,
        "Q": round(dl.Q, 3), "df": dl.df, "p_Q": round(dl.p_Q, 4),
        "I2": round(dl.I2, 1), "tau2": round(pm["tau2"], 5),
        "hk": {"df": pm["df"], "t_crit": round(pm["t_crit"], 3),
               "se": round(pm["se"], 5), "se_classic": round(pm["se_classic"], 5),
               "floor_applied": pm["floor_applied"]},
        "sens_dl": {"hr": round(math.exp(dl.est_random), 3),
                    "ci": _exp_ci(dl.ci_random), "tau2": round(dl.tau2, 5)},
        "sens_bayes": {"hr": round(math.exp(bay["est"]), 3), "ci": _exp_ci(bay["ci"]),
                       "tau_median": round(bay["tau_median"], 4),
                       "prior": f"half-normal(0, {bay['prior_tau']}) on tau"},
        "studies": per_study, "leave_one_out": loo,
        "total_n": sum(studies[s]["n_rand_exp"] + studies[s]["n_rand_ctrl"] for s in names),
    }


def risk_ratio(ev1, n1, ev0, n0):
    """Risk ratio with Katz log SE (no zero cells in these data)."""
    if min(ev1, ev0) == 0:
        return {}
    rr = (ev1 / n1) / (ev0 / n0)
    se = math.sqrt(1 / ev1 - 1 / n1 + 1 / ev0 - 1 / n0)
    lo, hi = math.exp(math.log(rr) - 1.959964 * se), math.exp(math.log(rr) + 1.959964 * se)
    return {"rr": round(rr, 3), "lo": round(lo, 3), "hi": round(hi, 3),
            "se_log": round(se, 4), "ev_exp": ev1, "n_exp": n1,
            "ev_ctrl": ev0, "n_ctrl": n0}


def binary_block(binary, studies):
    out = {}
    for key, o in binary.items():
        rows = [{"study": name, "label": studies[name]["label"],
                 **risk_ratio(d["ev_exp"], d["n_exp"], d["ev_ctrl"], d["n_ctrl"]),
                 "src": d["src"]}
                for name, d in o["data"].items()]
        out[key] = {"label": o["label"], "k": len(rows), "studies": rows,
                    "poolable": len(rows) >= 2}
    return out


def presentation(key, oc):
    """How each outcome is presented in the paper."""
    if oc["k"] < 3:
        mode = "trial-specific estimates, not pooled (two trials)"
    elif key in ("efficacy_primary", "biochemical_failure"):
        mode = "trial-specific estimates; pooled value is a descriptive summary"
    else:
        mode = "pooled effect"
    return {"presentation": mode, "k": oc["k"]}


def late_tissue_dose(studies, organ_dose, hr_outcomes, ab_values=(3.0, 5.0)):
    out = []
    for name, s in studies.items():
        od = organ_dose[name]
        gu = hr_outcomes["late_gu_g2"]["data"].get(name, {})
        gi = hr_outcomes["late_gi_g2"]["data"].get(name, {})
        for ab in ab_values:
            ctrl = eqd2(s["dose_ctrl_gy"], s["nfx_ctrl"], ab)
            row = {"study": name, "label": s["label"], "alpha_beta": ab,
                   "eqd2_control": round(ctrl, 1), "src_dose": od["src"]}
            for organ in ("urethra", "rectum"):
                e = eqd2(od[organ], s["nfx_exp"], ab)
                row[f"eqd2_{organ}"] = round(e, 1)
                row[f"delta_{organ}"] = round(e - ctrl, 1)
            row.update({"gu_hr": gu.get("hr"), "gu_lo": gu.get("lo"), "gu_hi": gu.get("hi"),
                        "gi_hr": gi.get("hr"), "gi_lo": gi.get("lo"), "gi_hi": gi.get("hi")})
            out.append(row)
    return out


def absolute_effects(risks, studies):
    out = []
    for label, study, horizon, pe, pc, src in risks:
        rd = (pe - pc) / 100.0
        n = 1.0 / abs(rd) if rd != 0 else None
        out.append({
            "outcome": label, "study": study, "label": studies[study]["label"],
            "horizon_y": horizon, "risk_exp_pct": pe, "risk_ctrl_pct": pc,
            "risk_difference_pp": round(rd * 100, 1),
            "measure": "NNH" if rd > 0 else ("NNT" if rd < 0 else "none"),
            "n": round(n) if n else None,
            "n_exact": round(n, 1) if n else None,
            "favours": "comparator" if rd > 0 else ("experimental" if rd < 0 else "neither"),
            "src": src})
    return out


def dose_effect(studies, hr_outcomes):
    out = []
    for name, s in studies.items():
        top = s.get("internal_boost") or s.get("dose_cap_gy") or s["dose_exp_gy"]
        eff = hr_outcomes["efficacy_primary"]["data"].get(name, {})
        gu = hr_outcomes["late_gu_g2"]["data"].get(name, {})
        out.append({
            "study": name, "label": s["label"], "top_dose_gy": top, "n_fx": s["nfx_exp"],
            "eqd2_top_ab1_5": round(eqd2(top, s["nfx_exp"], 1.5), 1),
            "dose_definition": ("40 Gy to 95% of CTV (integrated boost)"
                                if s.get("internal_boost") else
                                ("38.78 Gy max PTV (capped)" if s.get("dose_cap_gy")
                                 else "prescription dose")),
            "efficacy_hr": eff.get("hr"), "efficacy_lo": eff.get("lo"),
            "efficacy_hi": eff.get("hi"),
            "gu_hr": gu.get("hr"), "gu_lo": gu.get("lo"), "gu_hi": gu.get("hi")})
    return out


def information_size(oc, outcome_data, hr_target=0.85, alpha=0.05, power=0.80):
    """Required events (Schoenfeld, 1:1) and its diversity (D2) adjustment."""
    z_a = float(stats.norm.ppf(1 - alpha / 2.0))
    z_b = float(stats.norm.ppf(power))
    req = 4.0 * (z_a + z_b) ** 2 / (math.log(hr_target) ** 2)
    yi = np.array([s["yi"] for s in oc["studies"]])
    sei = np.array([s["sei"] for s in oc["studies"]])
    r = pooling.dersimonian_laird(yi, sei)
    v_f, v_r = r.se_fixed ** 2, r.se_random ** 2
    d2 = max(0.0, (v_r - v_f) / v_r) if v_r > 0 else 0.0
    adj = req / (1.0 - d2) if d2 < 1 else float("inf")

    accrued = []
    for st in oc["studies"]:
        d = outcome_data[st["study"]]
        e1, e0 = d.get("ev_exp"), d.get("ev_ctrl")
        accrued.append({"study": st["study"], "label": st["label"],
                        "events": e1 + e0 if e1 is not None and e0 is not None else None,
                        "src": d["src"]})
    known = [a["events"] for a in accrued if a["events"] is not None]
    total = sum(known)
    return {
        "hr_target": hr_target, "alpha": alpha, "power": power,
        "required_events": round(req),
        "diversity_D2_pct": round(d2 * 100, 1),
        "required_events_adjusted": round(adj) if math.isfinite(adj) else None,
        "accrued": accrued, "accrued_events_total": total,
        "n_studies_with_events": len(known), "k": oc["k"],
        "fraction_naive": round(total / req * 100, 1),
        "fraction_adjusted": (round(total / adj * 100, 1)
                              if math.isfinite(adj) and adj > 0 else None),
    }


def run(trials):
    S, HR = trials.STUDIES, trials.HR_OUTCOMES
    outcomes = {k: pool_outcome(k, HR[k], S) for k in POOLED_OUTCOMES}
    # Two-trial outcomes are presented trial by trial: no pooled estimate is kept,
    # only the heterogeneity statistics shown in Figure 4.
    for o in outcomes.values():
        if o["k"] < 3:
            for key in [x for x in o if x.startswith(("hr_", "ci_", "pred_interval", "sens_", "hk",
                                                          "leave_one_out", "estimator"))]:
                del o[key]
    return {
        "eqd2": eqd2_table(S),
        "outcomes": outcomes,
        "binary": binary_block(trials.BINARY_OUTCOMES, S),
        "presentation": {k: presentation(k, o) for k, o in outcomes.items()},
        "late_tissue_dose": late_tissue_dose(S, trials.ORGAN_DOSE, HR),
        "absolute_effects": absolute_effects(trials.ABSOLUTE_RISKS, S),
        "dose_effect": dose_effect(S, HR),
        "information_size": information_size(outcomes["efficacy_primary"],
                                             HR["efficacy_primary"]["data"]),
    }
