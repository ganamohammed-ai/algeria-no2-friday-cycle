#!/usr/bin/env python3
"""Check of the Friday intervals: Newey-West on observations (lag 7, as in friday_cycle.py) versus the calendar-day grid
(bandwidth 7, 28, 90 days). Also a placebo: the same model with a dummy for Mon/Tue/Wed/Thu (one at a time, Fri and Sat kept in the model).
Model: weather + regional background (adjust=2), as the ring-controlled values of the paper.
Usage: python detect_friday_se_check.py OUTDIR study_table.csv [cities_table.csv]   (study table includes 2019-2025; Ramadan 2026 appended)"""
import sys, os, math, csv
import numpy as np
import friday_cycle as F
A = F.A; D = F.D
F.RAM.append((D("2026-02-18"), D("2026-03-19")))

def nw_cal(X, res, dates, bw):
    d0 = dates[0]; idx = np.array([(d - d0).days for d in dates])
    G = np.zeros((idx.max() + 1, X.shape[1])); Xe = X * res[:, None]; G[idx] = Xe
    S = Xe.T @ Xe
    for l in range(1, bw + 1):
        Gl = G[l:].T @ G[:-l]; S += (1 - l / (bw + 1.0)) * (Gl + Gl.T)
    Xi = np.linalg.inv(X.T @ X)
    return Xi @ S @ Xi

def run(rows, site, placebo=None):
    sub = [r for r in rows if r["site"] == site]
    y, X, names, n = F.build(sub, 2, "fs")
    keys = ["no2_inner", "no2_ring", "blh", "wspd10", "t2m"]
    rr = [r for r in sub if r["usable"] == 1 and all(np.isfinite([r[k] for k in keys])) and r["blh"] > 0 and r["wspd10"] > 0 and r["no2_inner"] > 0 and r["no2_ring"] > 0]
    rr.sort(key=lambda r: r["date"]); dd = [D(r["date"]) for r in rr]
    if placebo is not None:
        X = np.column_stack([X, np.array([1.0 if d.weekday() == placebo else 0.0 for d in dd])]); names = names + ["placebo"]
    keep = [i for i in range(X.shape[1]) if i == 0 or X[:, i].std() > 0]
    X = X[:, keep]; names = [names[i] for i in keep]
    b, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ b
    target = "placebo" if placebo is not None else "Fri"
    i = names.index(target)
    out = {}
    V = {"nw7rows": A.newey_west(X, res, A.NW_LAG), "cal7": nw_cal(X, res, dd, 7), "cal28": nw_cal(X, res, dd, 28), "cal90": nw_cal(X, res, dd, 90)}
    for k, v in V.items():
        s = math.sqrt(v[i, i]); out[k] = (F.pct(b[i]), F.pct(b[i] - 1.96 * s), F.pct(b[i] + 1.96 * s), s)
    return out

def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True); rows = []
    for t in sys.argv[2:]: rows += A.read_csv(t)
    sites = [s for s in dict.fromkeys(r["site"] for r in rows) if s in F.STUDY or s in ("annaba", "tlemcen", "biskra", "djelfa", "chlef", "tebessa", "tiaret")]
    fri, plc = [], []
    for s in sites:
        o = run(rows, s)
        for k, (e, lo, hi, se) in o.items():
            fri.append(dict(site=s, se_method=k, est_pct=e, lo_pct=lo, hi_pct=hi, excludes_zero=int(hi < 0 or lo > 0), se_log=se))
        for wd in (0, 1, 2, 3):
            o = run(rows, s, placebo=wd)
            for k, (e, lo, hi, se) in o.items():
                plc.append(dict(site=s, weekday=wd, se_method=k, est_pct=e, excludes_zero=int(hi < 0 or lo > 0)))
    A.write_csv(os.path.join(out, "detect_friday_se_check.csv"), fri)
    A.write_csv(os.path.join(out, "detect_friday_placebo_days.csv"), plc)
main()
