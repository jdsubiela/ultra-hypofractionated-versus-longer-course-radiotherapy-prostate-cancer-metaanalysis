# Ultra-hypofractionated (stereotactic body) radiotherapy versus moderately hypofractionated or conventionally fractionated radiotherapy for localised prostate cancer: a systematic review and meta-analysis of randomised trials

This repository reproduces the numerical results of the review from the data extracted from the three trials that contribute to the synthesis (HYPO-RT-PC, PACE-B, NRG-GU005) and from the digitized Kaplan-Meier curves. The aggregate-data synthesis pools log hazard ratios with a Paule-Mandel estimate of between-trial variance and a modified Hartung-Knapp interval, with DerSimonian-Laird and a Bayesian normal-normal model as sensitivity analyses, and adds leave-one-out estimates, risk ratios for single-trial binary outcomes, equivalent doses in 2 Gy fractions, absolute effects and the required information size. Individual patient data are reconstructed from the published curves and numbers at risk with the Guyot (2012) algorithm, checked against the survival probabilities each trial prints, and analysed with Cox models at 3 and 5 years and with restricted mean survival time differences, pooled across trials. The analysis is deterministic; the seed (20260819) is fixed for completeness.

## How to run

```
pip install -r requirements.txt
python run_all.py
```

Results are written to `results/`:

- `synthesis.json`: aggregate-data meta-analysis
- `reconstruction_qa.json`: validation of each reconstructed arm
- `synthesis_ipd.json`: analyses on the reconstructed individual patient data
- `reconstructed_ipd/`: reconstructed individual patient data, regenerated on every run (not tracked)

## Layout

```
run_all.py              single entry point
data/trials.py          extracted trial results, each with its source
data/km_curves.py       numbers at risk, event totals and printed survival probabilities
data/digitized/         digitized Kaplan-Meier coordinates (time in years, survival in %)
sbrtma/pooling.py       Paule-Mandel, Hartung-Knapp, DerSimonian-Laird, Bayesian pooling
sbrtma/aggregate.py     aggregate-data synthesis
sbrtma/guyot.py         Guyot algorithm for IPD reconstruction
sbrtma/reconstruction.py  reconstruction of every arm and validation checks
sbrtma/ipd_analysis.py  Cox, landmark and RMST analyses on the reconstructed data
```

## Software

Python 3.14.6, numpy 2.4.2, scipy 1.17.1, pandas 2.3.3, lifelines 0.30.1.
