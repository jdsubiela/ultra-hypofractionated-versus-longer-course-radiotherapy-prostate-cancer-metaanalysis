"""Published Kaplan-Meier inputs for IPD reconstruction.

For each arm: the numbers-at-risk table (time in years), the total number of
events (used only when the curve is not truncated), the truncation time t_max,
and the survival probabilities (%) the trial prints, which serve as an
independent validation check (tolerance 0.5 percentage points).

Digitized coordinates are in data/digitized/<curve_id>.csv (time in years,
survival in %), extracted from the vector paths of the published figures.
"""

HYPO = "Nilsson 2026, Lancet Oncol"
PACEB = "van As 2024, N Engl J Med"
GU005 = "Ellis 2026, JAMA"

YEARS_0_15 = list(range(0, 16))
HALF_YEARS_0_5 = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

CURVES = {
    "hypo_ffs_uhf": {
        "trial": "HYPO-RT-PC", "arm": "experimental", "outcome": "Failure-free survival",
        "maxy": 100.0, "t_max": 10.0,
        "t_risk": YEARS_0_15,
        "n_risk": [589, 569, 543, 511, 489, 464, 432, 410, 377, 348, 266, 195,
                   140, 87, 60, 35],
        "reported_5y": 84.6, "reported_10y": 71.5,
        "src": f"{HYPO}, Figure 1 (numbers at risk); Suppl. Table 3 (5-y, 10-y)",
    },
    "hypo_ffs_cf": {
        "trial": "HYPO-RT-PC", "arm": "control", "outcome": "Failure-free survival",
        "maxy": 100.0, "t_max": 10.0,
        "t_risk": YEARS_0_15,
        "n_risk": [591, 580, 558, 532, 501, 474, 439, 404, 350, 309, 230, 171,
                   125, 88, 57, 39],
        "reported_5y": 85.5, "reported_10y": 65.0,
        "src": f"{HYPO}, Figure 1 (numbers at risk); Suppl. Table 3 (5-y, 10-y)",
    },
    "hypo_os_uhf": {
        "trial": "HYPO-RT-PC", "arm": "experimental", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": 10.0,
        "t_risk": YEARS_0_15,
        "n_risk": [589, 589, 581, 575, 566, 556, 541, 528, 510, 495, 405, 317,
                   238, 166, 114, 72],
        "reported_5y": 94.4, "reported_10y": 81.0,
        "src": f"{HYPO}, Figure 2 (numbers at risk); Suppl. Table 3 (5-y, 10-y)",
    },
    "hypo_os_cf": {
        "trial": "HYPO-RT-PC", "arm": "control", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": 10.0,
        "t_risk": YEARS_0_15,
        "n_risk": [591, 589, 585, 581, 572, 565, 546, 531, 508, 483, 395, 308,
                   234, 173, 118, 77],
        "reported_5y": 95.6, "reported_10y": 78.9,
        "src": f"{HYPO}, Figure 2 (numbers at risk); Suppl. Table 3 (5-y, 10-y)",
    },
    "paceb_ffbcf_sbrt": {
        "trial": "PACE-B", "arm": "experimental",
        "outcome": "Freedom from biochemical or clinical failure",
        "maxy": 100.0, "t_max": 6.0,
        "t_risk": [0, 1, 2, 3, 4, 5, 6],
        "n_risk": [433, 418, 405, 396, 380, 352, 217],
        "reported_5y": 95.8,
        "src": f"{PACEB}, Figure 1A (numbers at risk); Abstract (5-y)",
    },
    "paceb_ffbcf_ctrl": {
        "trial": "PACE-B", "arm": "control",
        "outcome": "Freedom from biochemical or clinical failure",
        "maxy": 100.0, "t_max": 6.0,
        "t_risk": [0, 1, 2, 3, 4, 5, 6],
        "n_risk": [441, 422, 411, 403, 387, 352, 206],
        "reported_5y": 94.6,
        "src": f"{PACEB}, Figure 1A (numbers at risk); Abstract (5-y)",
    },
    "paceb_os_sbrt": {
        "trial": "PACE-B", "arm": "experimental", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": 6.0,
        "t_risk": [0, 1, 2, 3, 4, 5, 6],
        "n_risk": [433, 426, 417, 408, 399, 385, 261],
        "src": f"{PACEB}, Figure 1D (numbers at risk)",
    },
    "paceb_os_ctrl": {
        "trial": "PACE-B", "arm": "control", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": 6.0,
        "t_risk": [0, 1, 2, 3, 4, 5, 6],
        "n_risk": [441, 425, 421, 417, 408, 396, 256],
        "src": f"{PACEB}, Figure 1D (numbers at risk)",
    },
    "gu005_dfs_sbrt": {
        "trial": "NRG-GU005", "arm": "experimental", "outcome": "Disease-free survival",
        "maxy": 100.0, "t_max": None,
        "t_risk": HALF_YEARS_0_5,
        "n_risk": [353, 349, 344, 325, 313, 268, 211, 152, 120, 88, 60],
        "events": 54,
        "reported_3y": 88.6,
        "src": f"{GU005}, Figure 2A (numbers at risk); Figure 2 legend (events); "
               "Results (3-y)",
    },
    "gu005_dfs_mhimrt": {
        "trial": "NRG-GU005", "arm": "control", "outcome": "Disease-free survival",
        "maxy": 100.0, "t_max": None,
        "t_risk": HALF_YEARS_0_5,
        "n_risk": [345, 325, 323, 314, 303, 255, 202, 139, 103, 80, 47],
        "events": 36,
        "reported_3y": 92.1,
        "src": f"{GU005}, Figure 2A (numbers at risk); Figure 2 legend (events); "
               "Results (3-y)",
    },
    "gu005_os_sbrt": {
        "trial": "NRG-GU005", "arm": "experimental", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": None,
        "t_risk": HALF_YEARS_0_5,
        "n_risk": [353, 349, 349, 343, 338, 292, 227, 166, 133, 97, 68],
        "events": 16,
        "reported_5y": 91.1,
        "src": f"{GU005}, Figure 2B (numbers at risk); Figure 2 legend (deaths); "
               "Table 3 (5-y)",
    },
    "gu005_os_mhimrt": {
        "trial": "NRG-GU005", "arm": "control", "outcome": "Overall survival",
        "maxy": 100.0, "t_max": None,
        "t_risk": HALF_YEARS_0_5,
        "n_risk": [345, 325, 323, 318, 309, 264, 211, 149, 113, 89, 55],
        "events": 13,
        "reported_5y": 94.5,
        "src": f"{GU005}, Figure 2B (numbers at risk); Figure 2 legend (deaths); "
               "Table 3 (5-y)",
    },
}

# Experimental and control arm of each reconstructed comparison.
PAIRS = [
    ("HYPO-RT-PC", "Failure-free survival", "hypo_ffs_uhf", "hypo_ffs_cf"),
    ("HYPO-RT-PC", "Overall survival", "hypo_os_uhf", "hypo_os_cf"),
    ("PACE-B", "Freedom from biochemical/clinical failure",
     "paceb_ffbcf_sbrt", "paceb_ffbcf_ctrl"),
    ("PACE-B", "Overall survival", "paceb_os_sbrt", "paceb_os_ctrl"),
    ("NRG-GU005", "Disease-free survival", "gu005_dfs_sbrt", "gu005_dfs_mhimrt"),
    ("NRG-GU005", "Overall survival", "gu005_os_sbrt", "gu005_os_mhimrt"),
]
