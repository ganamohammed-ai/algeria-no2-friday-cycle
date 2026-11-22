#!/usr/bin/env python3
"""Analysis B of analysis_plan_v2_detectability.md: 23 Mar-13 Jun window in 2020..2025, fitted within Jan 1-Aug 31 of each year."""
import datetime as dt, math, os, sys
import numpy as np
import lockdown_normalisation as L
import analyze_no2_ventilation as A
D = L.D
RAM = {2019:("05-06","06-04"),2020:("04-24","05-23"),2021:("04-13","05-12"),2022:("04-02","05-01"),
       2023:("03-23","04-20"),2024:("03-11","04-09"),2025:("03-01","03-29")}
STEPS = L.VERSION_STEPS

def fit_year(rows, Y, adjust):
    lo, hi = D(f"{Y}-01-01"), D(f"{Y}-08-31")
    rr = L.prep(rows, lo, hi)
    dates = [r["date"] for r in rr]
    y = np.log([r["no2_inner"] for r in rr])
    win = (D(f"{Y}-03-23"), D(f"{Y}-06-13"))
    x = L.in_window(dates, win)
    if x.sum() < 15: return None
    a, b = RAM[Y]; ra = (D(f"{Y}-{a}"), D(f"{Y}-{b}"))
    extra = []
    if adjust:
        extra += [np.log([r["wspd10"] for r in rr]), np.log([r["blh"] for r in rr]), np.array([r["t2m"] for r in rr])]
    extra.append(np.array([1.0 if ra[0] <= D(d) <= ra[1] else 0.0 for d in dates]))
    X, _ = A.design(dates, x, extra, with_month=True)
    X = X[:, [0, 1] + [j for j in range(2, X.shape[1]) if np.any(X[:, j] != 0)]]
    for s in STEPS:  # add a version step only if it adds rank (a step on 1 Jul equals the July dummy)
        if lo < s <= hi:
            c = np.array([1.0 if D(d) >= s else 0.0 for d in dates])
            X2 = np.column_stack([X, c])
            if np.linalg.matrix_rank(X2) > np.linalg.matrix_rank(X):
                X = X2
    keep = []
    for j in range(X.shape[1]):  # greedy: drop any column linearly dependent on earlier ones
        if np.linalg.matrix_rank(X[:, keep + [j]]) > len(keep):
            keep.append(j)
    X = X[:, keep]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    V = A.newey_west(X, y - X @ beta, A.NW_LAG)
    l, s = float(beta[1]), math.sqrt(V[1, 1]); p = lambda v: 100*(math.exp(v)-1)
    return dict(est=p(l), lo=p(l-1.96*s), hi=p(l+1.96*s), n_win=int(x.sum()), n=len(rr))


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rows = A.read_csv(args.table) + A.read_csv(args.controls)
    out = []
    for s in L.STUDY + L.CONTROLS:
        sr = [r for r in rows if r["site"] == s]
        for adj in (False, True):
            for Y in range(2020, 2026):
                f = fit_year(sr, Y, adj)
                if f: out.append(dict(site=s, model="adjusted" if adj else "unadjusted", year=Y, **f))
    A.write_csv(os.path.join(args.out, "detect_seasonal.csv"), out)
    print("site".ljust(17), *[str(y).rjust(8) for y in range(2020,2026)], " rank2020")
    for s in L.STUDY + L.CONTROLS:
        r = [o for o in out if o["site"]==s and o["model"]=="adjusted"]
        v = {o["year"]:o["est"] for o in r}
        if 2020 in v:
            rank = sum(1 for y in range(2021,2026) if y in v and abs(v[y]) >= abs(v[2020]))
            print(s.ljust(17), *[f"{v.get(y,float('nan')):+.1f}".rjust(8) for y in range(2020,2026)], f" {rank}/5 springs >= |2020|")
    print("\nunadjusted")
    for s in L.STUDY + L.CONTROLS:
        v = {o["year"]:o["est"] for o in out if o["site"]==s and o["model"]=="unadjusted"}
        print(s.ljust(17), *[f"{v.get(y,float('nan')):+.1f}".rjust(8) for y in range(2020,2026)])

if __name__ == "__main__":
    main()
