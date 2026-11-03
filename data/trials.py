"""Aggregate data extracted from the three included trials.

Hazard ratios are experimental (ultra-hypofractionated) versus control; HR < 1
favours the ultra-hypofractionated arm. Every record carries `src`, a short
citation of where the value is printed.

Sources:
  HYPO-RT-PC  Nilsson 2026, Lancet Oncol 27:293-301 (+ supplementary appendix)
  PACE-B      van As 2024, N Engl J Med 391:1413-25 (+ supplementary appendix)
  NRG-GU005   Ellis 2026, JAMA (+ Supplement 2)
"""

HYPO = "Nilsson 2026, Lancet Oncol"
PACEB = "van As 2024, N Engl J Med"
GU005 = "Ellis 2026, JAMA"

# Schedules and sample sizes. dose_ctrl_gy/nfx_ctrl for PACE-B is the
# predominant control schedule delivered (62 Gy / 20 fx).
STUDIES = {
    "HYPO-RT-PC": {
        "label": "HYPO-RT-PC (10-y)",
        "comparator_type": "conventional",
        "n_rand_exp": 598, "n_rand_ctrl": 602,
        "dose_exp_gy": 42.7, "nfx_exp": 7,
        "dose_ctrl_gy": 78.0, "nfx_ctrl": 39,
        "internal_boost": None,
        "dose_cap_gy": None,
        "src": f"{HYPO}, Abstract; Suppl. Table 1",
    },
    "PACE-B": {
        "label": "PACE-B (5-y)",
        "comparator_type": "mixed",
        "n_rand_exp": 433, "n_rand_ctrl": 441,
        "dose_exp_gy": 36.25, "nfx_exp": 5,
        "dose_ctrl_gy": 62.0, "nfx_ctrl": 20,
        "internal_boost": 40.0,   # secondary dose to 95% of the CTV
        "dose_cap_gy": None,
        "src": f"{PACEB}, Abstract; Methods; Table 1",
    },
    "NRG-GU005": {
        "label": "NRG-GU005 (3-y)",
        "comparator_type": "moderate",
        "n_rand_exp": 353, "n_rand_ctrl": 345,
        "dose_exp_gy": 36.25, "nfx_exp": 5,
        "dose_ctrl_gy": 70.0, "nfx_ctrl": 28,
        "internal_boost": None,
        "dose_cap_gy": 38.78,     # maximum PTV dose
        "src": f"{GU005}, Abstract; Methods; Table 1",
    },
}

HR_OUTCOMES = {
    "efficacy_primary": {
        "label": "Composite oncological failure (trial-defined primary)",
        "data": {
            "HYPO-RT-PC": {"hr": 0.85, "lo": 0.70, "hi": 1.04, "ci": 95,
                           "ev_exp": 178, "n_exp": 589, "ev_ctrl": 205, "n_ctrl": 591,
                           "src": f"{HYPO}, Suppl. Table 3 (unadjusted)"},
            "PACE-B": {"hr": 0.73, "lo": 0.44, "hi": 1.21, "ci": 95,
                       "ev_exp": 26, "n_exp": 433, "ev_ctrl": 36, "n_ctrl": 441,
                       "src": f"{PACEB}, Results (95% CI)"},
            "NRG-GU005": {"hr": 1.38, "lo": 0.91, "hi": 2.09, "ci": 95,
                          "ev_exp": 54, "n_exp": 353, "ev_ctrl": 36, "n_ctrl": 345,
                          "src": f"{GU005}, Table 3 (unadjusted); events, Figure 2 legend"},
        },
    },
    "biochemical_failure": {
        "label": "Biochemical failure",
        "data": {
            "HYPO-RT-PC": {"hr": 0.86, "lo": 0.70, "hi": 1.06, "ci": 95,
                           "ev_exp": 174, "n_exp": 589, "ev_ctrl": 198, "n_ctrl": 591,
                           "src": f"{HYPO}, Suppl. Table 3 (unadjusted)"},
            "PACE-B": {"hr": 0.73, "lo": 0.44, "hi": 1.21, "ci": 95,
                       "ev_exp": 26, "n_exp": 433, "ev_ctrl": 36, "n_ctrl": 441,
                       "src": f"{PACEB}, Results (biochemical or clinical failure)"},
            "NRG-GU005": {"hr": 1.82, "lo": 1.03, "hi": 3.33, "ci": 95,
                          "ev_exp": None, "n_exp": 353, "ev_ctrl": None, "n_ctrl": 345,
                          "src": f"{GU005}, Table 3 (unadjusted)"},
        },
    },
    "late_gu_g2": {
        "label": "Late genitourinary adverse events, grade >=2",
        "data": {
            "HYPO-RT-PC": {"hr": 1.01, "lo": 0.81, "hi": 1.25, "ci": 95,
                           "src": f"{HYPO}, Abstract"},
            "PACE-B": {"hr": 1.59, "lo": 1.18, "hi": 2.12, "ci": 95,
                       "src": f"{PACEB}, Results"},
        },
    },
    "late_gi_g2": {
        "label": "Late gastrointestinal adverse events, grade >=2",
        "data": {
            "HYPO-RT-PC": {"hr": 0.94, "lo": 0.70, "hi": 1.28, "ci": 95,
                           "src": f"{HYPO}, Abstract"},
            "PACE-B": {"hr": 1.03, "lo": 0.68, "hi": 1.56, "ci": 95,
                       "src": f"{PACEB}, Results"},
        },
    },
    "overall_survival": {
        "label": "Overall survival",
        "data": {
            "PACE-B": {"hr": 1.41, "lo": 0.90, "hi": 2.20, "ci": 95,
                       "ev_exp": 46, "n_exp": 433, "ev_ctrl": 33, "n_ctrl": 441,
                       "src": f"{PACEB}, Results"},
            "HYPO-RT-PC": {"hr": 0.91, "lo": 0.75, "hi": 1.11, "ci": 95,
                           "ev_exp": 189, "n_exp": 589, "ev_ctrl": 208, "n_ctrl": 591,
                           "src": f"{HYPO}, Suppl. Table 3 (unadjusted)"},
            "NRG-GU005": {"hr": 1.21, "lo": 0.54, "hi": 2.33, "ci": 95,
                          "ev_exp": 16, "n_exp": 353, "ev_ctrl": 13, "n_ctrl": 345,
                          "src": f"{GU005}, Table 3 (unadjusted); deaths, Figure 2 legend"},
        },
    },
}

# Counts by arm (single-trial outcomes; risk ratios are computed downstream).
BINARY_OUTCOMES = {
    "gu_g3plus": {
        "label": "Grade >=3 genitourinary adverse events",
        "data": {
            "NRG-GU005": {"ev_exp": 2, "n_exp": 352, "ev_ctrl": 8, "n_ctrl": 323,
                          "src": f"{GU005}, Supplement 2, eTable 4"},
        },
    },
    "gi_g3plus": {
        "label": "Grade >=3 gastrointestinal adverse events",
        "data": {
            "NRG-GU005": {"ev_exp": 3, "n_exp": 352, "ev_ctrl": 7, "n_ctrl": 323,
                          "src": f"{GU005}, Supplement 2, eTable 4"},
        },
    },
    "epic_bowel_mcid_2y": {
        "label": "EPIC-26 bowel domain, MCID decline at 2 y",
        "data": {
            "NRG-GU005": {"ev_exp": 101, "n_exp": 289, "ev_ctrl": 113, "n_ctrl": 258,
                          "src": f"{GU005}, Table 2 (denominators back-calculated "
                                 "from printed percentages)"},
        },
    },
    "epic_urinary_irritative_mcid_2y": {
        "label": "EPIC-26 urinary irritative/obstructive, MCID decline at 2 y",
        "data": {
            "NRG-GU005": {"ev_exp": 102, "n_exp": 288, "ev_ctrl": 102, "n_ctrl": 303,
                          "src": f"{GU005}, Table 2 (denominators back-calculated "
                                 "from printed percentages)"},
        },
    },
}

# Prescription-level surrogate doses (Gy) for urethra (inside the CTV) and
# rectum (at the PTV edge) in the experimental arm.
ORGAN_DOSE = {
    "HYPO-RT-PC": {"urethra": 42.7, "rectum": 42.7, "src": f"{HYPO}, Suppl. Table 1"},
    "PACE-B": {"urethra": 40.0, "rectum": 36.25, "src": f"{PACEB}, Methods"},
    "NRG-GU005": {"urethra": 38.78, "rectum": 36.25, "src": f"{GU005}, Methods"},
}

# Risks by arm (%) at a stated horizon, for absolute effects.
# (outcome, study, horizon_y, risk_exp_pct, risk_ctrl_pct, src)
ABSOLUTE_RISKS = [
    ("Late GU adverse events, grade >=2", "HYPO-RT-PC", 10, 28.3, 30.1,
     f"{HYPO}, Suppl. Table 4"),
    ("Late GU adverse events, grade >=2", "PACE-B", 5, 26.9, 18.3,
     f"{PACEB}, Results"),
    ("Late GI adverse events, grade >=2", "HYPO-RT-PC", 10, 14.4, 14.5,
     f"{HYPO}, Suppl. Table 4"),
    ("Late GI adverse events, grade >=2", "PACE-B", 5, 10.7, 10.2,
     f"{PACEB}, Results"),
    ("GU adverse events, grade >=3", "NRG-GU005", 3, 0.57, 2.48,
     f"{GU005}, Supplement 2, eTable 4 (2/352 vs 8/323)"),
    ("EPIC-26 bowel domain, MCID decline", "NRG-GU005", 2, 34.9, 43.8,
     f"{GU005}, Table 2"),
    ("Oncological failure", "HYPO-RT-PC", 10, 28.5, 35.0,
     f"{HYPO}, Suppl. Table 3 (100 minus 10-y FFS)"),
    ("Oncological failure", "PACE-B", 5, 4.2, 5.4,
     f"{PACEB}, Abstract (100 minus 5-y freedom from failure)"),
    ("Oncological failure", "NRG-GU005", 3, 11.4, 7.9,
     f"{GU005}, Results (100 minus 3-y DFS)"),
]
