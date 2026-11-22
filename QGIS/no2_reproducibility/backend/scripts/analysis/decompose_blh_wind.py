#!/usr/bin/env python3
"""Split the ventilation product into its two parts.

VC10 = BLH x |wind at 10 m|.  A well-mixed box model predicts that a vertically
integrated column depends on the wind (advection) but not directly on BLH.  This
script fits, for each site,

    ln(NO2 inner) ~ ln(BLH) + ln(wind10) + T2m + trend + month + weekday

with the same sample, controls and Newey-West errors (lag 7) as the main model
(analyze_no2_ventilation.py), and writes beta_BLH and beta_wind with their 95% CI.
It also reports the correlation of ln(BLH) with ln(wind10).

Usage:  python decompose_blh_wind.py --table no2_ventilation_daily_v2.csv --out outputs
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze_no2_ventilation as A  # noqa: E402

SITES = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"]


def fit_two(rows):
    keys = ["no2_inner", "blh", "wspd10", "t2m"]
    rr = [r for r in rows if r["usable"] == 1 and all(np.isfinite([r[k] for k in keys]))
          and r["blh"] > 0 and r["wspd10"] > 0 and r["no2_inner"] > 0]
    rr.sort(key=lambda r: r["date"])
    y = np.log([r["no2_inner"] for r in rr])
    lb = np.log([r["blh"] for r in rr])
    lw = np.log([r["wspd10"] for r in rr])
    t2 = np.array([r["t2m"] for r in rr])
    X, _ = A.design([r["date"] for r in rr], lb, [lw, t2], with_month=True)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    V = A.newey_west(X, res, A.NW_LAG)
    out = {"n": len(rr), "corr_lnBLH_lnWind": float(np.corrcoef(lb, lw)[0, 1])}
    for name, idx in (("blh", 1), ("wind10", 2)):
        se = math.sqrt(V[idx, idx])
        X0 = np.delete(X, idx, axis=1)
        b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
        r0 = y - X0 @ b0
        out["beta_" + name] = float(beta[idx])
        out["ci_low_" + name] = float(beta[idx] - 1.96 * se)
        out["ci_high_" + name] = float(beta[idx] + 1.96 * se)
        out["partial_R2_" + name] = float(1.0 - (res @ res) / (r0 @ r0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    rows = A.read_csv(a.table)
    results = []
    for site in SITES:
        sub = [r for r in rows if r["site"] == site]
        if not sub:
            continue
        f = fit_two(sub)
        f = dict(site=site, **f)
        results.append(f)
        print("%-15s n=%4d  beta_BLH=%6.3f (%6.3f..%6.3f) pR2=%.3f | beta_wind=%6.3f (%6.3f..%6.3f) pR2=%.3f | r(lnBLH,lnW)=%5.2f"
              % (site, f["n"], f["beta_blh"], f["ci_low_blh"], f["ci_high_blh"], f["partial_R2_blh"],
                 f["beta_wind10"], f["ci_low_wind10"], f["ci_high_wind10"], f["partial_R2_wind10"],
                 f["corr_lnBLH_lnWind"]))
    A.write_csv(os.path.join(a.out, "sensitivity_blh_wind_split.csv"), results)


if __name__ == "__main__":
    main()
