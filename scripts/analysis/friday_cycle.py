#!/usr/bin/env python3
"""Weekly cycle of the TROPOMI NO2 column after removing weather (Friday/Saturday weekend in Algeria).

Model per site (usable days, 2019-08-07..2025-12-31), OLS with Newey-West errors (lag 7):
  ln NO2 = a + [b1 ln wind10 + b2 ln BLH + b3 T2m] + trend + month + Ramadan + containment(2020-03-23..07-28)
           + weekday dummies (Mon..Sat; Sunday is the reference)
'unadjusted' leaves out the bracket. Effects in % relative to Sunday = 100*(exp(b)-1).
Also reports Friday and Saturday relative to the mean of Sun-Thu (second model with Fri, Sat dummies only)
and the Ramadan effect. Ramadan dates are approximate (+/- 1 day).
Usage: python friday_cycle.py --table T.csv --controls C.csv --out outputs
"""
import argparse, datetime as dt, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_no2_ventilation as A
D = dt.date.fromisoformat
STUDY = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"]
CONTROLS = ["control_tell", "control_highland"]
RAMADAN = [("2019-05-06", "2019-06-04"), ("2020-04-24", "2020-05-23"), ("2021-04-13", "2021-05-12"),
           ("2022-04-02", "2022-05-01"), ("2023-03-23", "2023-04-20"), ("2024-03-11", "2024-04-09"),
           ("2025-03-01", "2025-03-29")]
RAM = [(D(a), D(b)) for a, b in RAMADAN]
COV = (D("2020-03-23"), D("2020-07-28"))
DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]  # python weekday 0..5; Sunday (6) = reference


def pct(b): return 100.0 * (math.exp(b) - 1.0)


def build(rows, adjust, mode):
    keys = ["no2_inner", "no2_ring", "blh", "wspd10", "t2m"]
    rr = [r for r in rows if r["usable"] == 1 and all(np.isfinite([r[k] for k in keys]))
          and r["blh"] > 0 and r["wspd10"] > 0 and r["no2_inner"] > 0 and r["no2_ring"] > 0]
    rr.sort(key=lambda r: r["date"])
    dd = [D(r["date"]) for r in rr]
    y = np.log([r["no2_inner"] for r in rr])
    cols, names = [np.ones(len(dd))], ["const"]
    if adjust:
        cols += [np.log([r["wspd10"] for r in rr]), np.log([r["blh"] for r in rr]), np.array([r["t2m"] for r in rr])]
        names += ["ln_wind", "ln_blh", "t2m"]
    if adjust == 2:
        cols.append(np.log([r["no2_ring"] for r in rr])); names.append("ln_ring")
    cols.append(np.array([(d - D("2019-01-01")).days / 365.25 for d in dd])); names.append("trend")
    for m in range(2, 13):
        cols.append(np.array([1.0 if d.month == m else 0.0 for d in dd])); names.append("m%02d" % m)
    cols.append(np.array([1.0 if any(a <= d <= b for a, b in RAM) else 0.0 for d in dd])); names.append("ramadan")
    cols.append(np.array([1.0 if COV[0] <= d <= COV[1] else 0.0 for d in dd])); names.append("containment")
    if mode == "all":
        for k, nm in enumerate(DAYS):
            cols.append(np.array([1.0 if d.weekday() == k else 0.0 for d in dd])); names.append(nm)
    else:  # fri/sat vs Sun-Thu
        for k, nm in ((4, "Fri"), (5, "Sat")):
            cols.append(np.array([1.0 if d.weekday() == k else 0.0 for d in dd])); names.append(nm)
    return y, np.column_stack(cols), names, len(rr)


def fit(rows, adjust, mode):
    y, X, names, n = build(rows, adjust, mode)
    keep = [i for i in range(X.shape[1]) if i == 0 or X[:, i].std() > 0]  # drop all-zero dummies in short sub-samples
    X = X[:, keep]; names = [names[i] for i in keep]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    V = A.newey_west(X, res, A.NW_LAG)
    out = {}
    for nm in names:
        if nm in DAYS + ["Fri", "Sat", "ramadan"] and nm in names:
            i = names.index(nm); s = math.sqrt(V[i, i])
            out[nm] = (pct(b[i]), pct(b[i] - 1.96 * s), pct(b[i] + 1.96 * s))
    return n, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True); ap.add_argument("--controls", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    rows = A.read_csv(a.table) + A.read_csv(a.controls)
    res = []
    for s in STUDY + CONTROLS:
        sub = [r for r in rows if r["site"] == s]
        row = dict(site=s)
        for adj in (0, 1, 2):
            tag = {0: "raw", 1: "adj", 2: "ring"}[adj]
            n, o = fit(sub, adj, "fs")
            row["n"] = n
            for nm in ("Fri", "Sat", "ramadan"):
                row["%s_%s_pct" % (nm, tag)], row["%s_%s_lo" % (nm, tag)], row["%s_%s_hi" % (nm, tag)] = o[nm]
            n, o = fit(sub, adj, "all")
            for nm in DAYS:
                row["all_%s_%s_pct" % (nm, tag)] = o[nm][0]
        res.append(row)
        print("   ring-controlled Fri %6.1f (%6.1f..%6.1f)  Sat %6.1f (%6.1f..%6.1f)" % (row["Fri_ring_pct"], row["Fri_ring_lo"], row["Fri_ring_hi"], row["Sat_ring_pct"], row["Sat_ring_lo"], row["Sat_ring_hi"]))
        print("%-17s n=%4d  Fri raw %6.1f (%6.1f..%6.1f)  adj %6.1f (%6.1f..%6.1f) | Sat adj %6.1f (%6.1f..%6.1f) | Ramadan adj %6.1f (%6.1f..%6.1f)" % (
            s, row["n"], row["Fri_raw_pct"], row["Fri_raw_lo"], row["Fri_raw_hi"], row["Fri_adj_pct"], row["Fri_adj_lo"], row["Fri_adj_hi"],
            row["Sat_adj_pct"], row["Sat_adj_lo"], row["Sat_adj_hi"], row["ramadan_adj_pct"], row["ramadan_adj_lo"], row["ramadan_adj_hi"]))
    A.write_csv(os.path.join(a.out, "friday_cycle.csv"), res)


if __name__ == "__main__":
    main()
