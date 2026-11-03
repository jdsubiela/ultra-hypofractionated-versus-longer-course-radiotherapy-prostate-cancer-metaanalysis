#!/usr/bin/env python3
"""Reproduce all numerical results.

  python run_all.py

Writes:
  results/synthesis.json            aggregate-data meta-analysis
  results/reconstruction_qa.json    IPD reconstruction validation checks
  results/synthesis_ipd.json        analyses on the reconstructed IPD
  results/reconstructed_ipd/*.csv   reconstructed IPD (regenerated each run)
"""
import json
import os
import platform
import sys

import lifelines
import numpy as np
import pandas as pd
import scipy

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)  # also needed under python -I / -P

from data import km_curves, trials  # noqa: E402
from sbrtma import aggregate, ipd_analysis, reconstruction  # noqa: E402

SEED = 20260819
RESULTS = os.path.join(ROOT, "results")
IPD_DIR = os.path.join(RESULTS, "reconstructed_ipd")
DIGITIZED = os.path.join(ROOT, "data", "digitized")


def write(name, obj):
    with open(os.path.join(RESULTS, name), "w") as f:
        json.dump(obj, f, indent=2)


def r2(x):
    """Two decimals, rounding half up (1.055 is 1.06)."""
    from decimal import Decimal, ROUND_HALF_UP
    return str(Decimal(str(x)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def main():
    np.random.seed(SEED)
    os.makedirs(RESULTS, exist_ok=True)
    versions = {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__, "pandas": pd.__version__,
                "lifelines": lifelines.__version__}

    agg = aggregate.run(trials)
    write("synthesis.json", {"seed": SEED, "versions": versions, **agg})

    qa = reconstruction.reconstruct_all(km_curves.CURVES, DIGITIZED, IPD_DIR)
    write("reconstruction_qa.json", qa)

    ipd = ipd_analysis.run(km_curves.PAIRS, qa, IPD_DIR)
    write("synthesis_ipd.json", {"seed": SEED, "versions": versions, **ipd})

    print("Aggregate meta-analysis (Paule-Mandel + modified Hartung-Knapp)")
    for key, o in agg["outcomes"].items():
        est = (f"HR {r2(o['hr_random'])} ({r2(o['ci_random'][0])} to {r2(o['ci_random'][1])})"
               if "hr_random" in o else "not pooled")
        print(f"  {o['label']:<55} k={o['k']}  {est}  "
              f"I2={o['I2']:.0f}%  tau2={o['tau2']:.4f}  "
              f"[{agg['presentation'][key]['presentation']}]")
    n_pass = sum(m["pass"] for a in qa.values() for m in a["milestones"])
    n_all = sum(len(a["milestones"]) for a in qa.values())
    n_ok = sum(a["pass"] for a in qa.values())
    print(f"IPD reconstruction: {n_ok}/{len(qa)} arms pass validation; "
          f"{n_pass}/{n_all} milestone checks pass")
    print("IPD, common-effect HR with follow-up truncated at each horizon")
    for g, blocks in ipd["meta"].items():
        for key, b in blocks.items():
            if key.startswith("landmark"):
                print(f"  {g:<18} {key:<12} HR {b['hr']:.2f} "
                      f"({b['lo']:.2f} to {b['hi']:.2f})  I2={b['I2']:.0f}%")
    print(f"Results written to {os.path.relpath(RESULTS, os.getcwd())}")


if __name__ == "__main__":
    main()
