#!/usr/bin/env python3
"""Does the ventilation elasticity change when TROPOMI processor-version steps are added?

Main model (analyze_no2_ventilation.fit_site): ln NO2 ~ ln VC10 + T2m + trend + month + weekday.
Here the same model plus step dummies at the operational start dates of v2.2.0 ... v2.7.1
(2021-07-01, 2021-11-14, 2022-07-17, 2023-03-12, 2023-11-26, 2024-09-08).
Usage: python sensitivity_version_steps.py --table T.csv --controls C.csv --out outputs
"""
import argparse, datetime as dt, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_no2_ventilation as A
D = dt.date.fromisoformat
STEPS = [D(x) for x in ("2021-07-01", "2021-11-14", "2022-07-17", "2023-03-12", "2023-11-26", "2024-09-08")]
SITES = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud", "control_tell", "control_highland"]

def fit(rows, steps):
    rr = [r for r in rows if r["usable"] == 1 and all(np.isfinite([r["no2_inner"], r["vc10"], r["t2m"]]))
          and r["vc10"] > 0 and r["no2_inner"] > 0]
    rr.sort(key=lambda r: r["date"])
    y = np.log([r["no2_inner"] for r in rr]); x = np.log([r["vc10"] for r in rr])
    dates = [r["date"] for r in rr]
    extra = [np.array([r["t2m"] for r in rr])]
    if steps:
        extra += [np.array([1.0 if D(d) >= s else 0.0 for d in dates]) for s in STEPS]
    X, names = A.design(dates, x, extra, with_month=True)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    V = A.newey_west(X, res, A.NW_LAG)
    se = math.sqrt(V[1, 1])
    return len(rr), b[1], b[1] - 1.96 * se, b[1] + 1.96 * se

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True); ap.add_argument("--controls", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    rows = A.read_csv(a.table) + A.read_csv(a.controls)
    out = []
    for s in SITES:
        sub = [r for r in rows if r["site"] == s]
        n0, b0, l0, h0 = fit(sub, False); n1, b1, l1, h1 = fit(sub, True)
        out.append(dict(site=s, n=n0, beta_main=b0, lo_main=l0, hi_main=h0, beta_with_steps=b1, lo_steps=l1, hi_steps=h1))
        print("%-17s n=%4d  %.3f (%.3f..%.3f)  ->  %.3f (%.3f..%.3f)" % (s, n0, b0, l0, h0, b1, l1, h1))
    A.write_csv(os.path.join(a.out, "sensitivity_version_steps.csv"), out)

if __name__ == "__main__":
    main()
