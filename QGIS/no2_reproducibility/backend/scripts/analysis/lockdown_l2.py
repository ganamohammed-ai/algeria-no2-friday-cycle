#!/usr/bin/env python3
"""Run the pre-set lockdown normalisation (lockdown_normalisation.py, sample A) on the single-version L2 table,
and compare the L2 column with the Level-3 column on the same days.

Study sites and, when their clips are present in the table, the two background controls (negative-control rule).
Nothing in the model is changed; only the NO2 input differs.

Usage: python scripts\\analysis\\lockdown_l2.py --l2 outputs\\l2_v240_daily_q75.csv --l3 outputs\\no2_ventilation_daily_v2.csv --out outputs --tag q75
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_no2_ventilation as A
import lockdown_normalisation as LN


def compare(l2, l3):
    k3 = {(r["site"], r["date"]): r for r in l3 if r["usable"] == 1}
    out = []
    for s in LN.STUDY + LN.CONTROLS:
        x, y = [], []
        for r in l2:
            if r["site"] != s or r["usable"] != 1:
                continue
            b = k3.get((s, r["date"]))
            if b and r["no2_inner"] > 0 and b["no2_inner"] > 0:
                x.append(np.log(r["no2_inner"])); y.append(np.log(b["no2_inner"]))
        if len(x) > 20:
            x, y = np.array(x), np.array(y)
            out.append(dict(site=s, n_common_days=len(x), corr_ln=float(np.corrcoef(x, y)[0, 1]),
                            median_ratio_L2_over_L3=float(np.exp(np.median(x - y))),
                            mean_ratio_L2_over_L3=float(np.exp(np.mean(x - y)))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l2", required=True); ap.add_argument("--l3", required=True, nargs="+")
    ap.add_argument("--out", default="."); ap.add_argument("--tag", default="q75")
    a = ap.parse_args()
    l2 = A.read_csv(a.l2)
    l3 = [r for t in a.l3 for r in A.read_csv(t)]
    L, P, M = [], [], []
    for label, names in (("study_L2", LN.STUDY), ("control_L2", LN.CONTROLS)):
        sd = {s: [r for r in l2 if r["site"] == s] for s in names}
        sd = {k: v for k, v in sd.items() if v}
        if not sd:
            continue
        L += LN.lockdown_table(sd, label); P += LN.placebo_table(sd, label); M += LN.ml_table(sd, label)
    for name, rows in (("l2_vs_l3_%s.csv" % a.tag, compare(l2, l3)), ("l2_lockdown_effect_%s.csv" % a.tag, L),
                       ("l2_lockdown_placebo_%s.csv" % a.tag, P), ("l2_lockdown_ml_%s.csv" % a.tag, M)):
        A.write_csv(os.path.join(a.out, name), rows); print("wrote", name, len(rows))


if __name__ == "__main__":
    main()
