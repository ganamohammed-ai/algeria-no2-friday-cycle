#!/usr/bin/env python3
"""Calibration of the inference used in Analysis A (detect_injection.py) and Table 8.
Same model, windows and injection as detect_injection.py; only the variance estimator / decision rule changes:
  nw7rows : Newey-West, lag 7 OBSERVATIONS (original)
  cal28   : Newey-West on the calendar-day grid, bandwidth 28 days (missing days contribute zero)
  cal90   : same, bandwidth 90 days (longer than the 83-day window)
  emp     : no standard error; decision by the 2.5th percentile of the SAME site's window estimates at d=0 (one-sided)
Usage: python detect_calibrate.py OUTDIR TABLE1 [TABLE2 ...]   (tables with columns site,date,no2_inner,wspd10,blh,t2m,...)"""
import sys, os, math, datetime as dt
import numpy as np
import lockdown_normalisation as L
import analyze_no2_ventilation as A
import detect_injection as I
D = L.D
WIN = L.WINDOWS["primary_0323_0613"]
DELTAS = [0.0, -0.05, -0.10, -0.15, -0.20, -0.30]

def nw_cal(X, res, dates, bw):
    d0 = D(dates[0]); idx = np.array([(D(d) - d0).days for d in dates])
    n_days = idx.max() + 1
    Xe = X * res[:, None]
    G = np.zeros((n_days, X.shape[1])); G[idx] = Xe
    S = Xe.T @ Xe
    for l in range(1, bw + 1):
        w = 1.0 - l / (bw + 1.0)
        Gl = G[l:].T @ G[:-l]
        S += w * (Gl + Gl.T)
    XtX_inv = np.linalg.inv(X.T @ X)
    return XtX_inv @ S @ XtX_inv

def fitv(rr, win, adjust, placebo_win=None):
    dates = [r["date"] for r in rr]
    y = np.log([r["no2_inner"] for r in rr])
    lb = np.log([r["blh"] for r in rr]); lw = np.log([r["wspd10"] for r in rr]); t2 = np.array([r["t2m"] for r in rr])
    real = L.in_window(dates, win)
    target = L.in_window(dates, placebo_win) if placebo_win is not None else real
    extra = ([lw, lb, t2] if adjust else []) + [L.ramadan(dates)]
    if placebo_win is not None: extra.append(real)
    X, _ = A.design(dates, target, extra, with_month=True)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    se = {"nw7rows": math.sqrt(A.newey_west(X, res, A.NW_LAG)[1, 1]),
          "cal28": math.sqrt(nw_cal(X, res, dates, 28)[1, 1]),
          "cal90": math.sqrt(nw_cal(X, res, dates, 90)[1, 1])}
    return float(beta[1]), se

pc = lambda x: 100 * (np.exp(x) - 1)
def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True); tabs = sys.argv[2:]
    rows = []
    for t in tabs: rows += A.read_csv(t)
    sites = sorted({r["site"] for r in rows}, key=lambda s: (s not in L.STUDY, s))
    resid, rows_out, true2020, nullrows = [], [], [], []
    for site in sites:
        rr = L.prep([r for r in rows if r["site"] == site], *L.SAMPLE_A)
        if len(rr) < 200: continue
        ws = I.windows(rr)
        for adjust in (False, True):
            mod = "adjusted" if adjust else "unadjusted"
            # 2020 estimate with each se
            l, se = fitv(rr, WIN, adjust)
            for v, s_ in se.items():
                true2020.append(dict(site=site, model=mod, se_method=v, est_pct=pc(l), lo_pct=pc(l - 1.96 * s_), hi_pct=pc(l + 1.96 * s_)))
            store = {}
            for d in DELTAS:
                E, S = [], {k: [] for k in ("nw7rows", "cal28", "cal90")}
                for w in ws:
                    try: l, se = fitv(I.inject(rr, w, d), WIN, adjust, placebo_win=w)
                    except np.linalg.LinAlgError: continue
                    E.append(l)
                    for k in S: S[k].append(se[k])
                store[d] = (np.array(E), {k: np.array(v) for k, v in S.items()})
            e0 = store[0.0][0]
            for w_, v_ in zip(ws, e0): nullrows.append(dict(site=site, model=mod, start=str(w_[0]), est_pct=float(pc(v_))))
            thr = np.percentile(e0, 2.5)
            thr2 = np.percentile(e0, [2.5, 97.5])
            for d in DELTAS:
                E, S = store[d]
                for v in ("nw7rows", "cal28", "cal90"):
                    lo, hi = E - 1.96 * S[v], E + 1.96 * S[v]
                    rows_out.append(dict(site=site, model=mod, method=v, d_true_pct=100 * d, n_windows=len(E),
                                         mean_est_pct=pc(E).mean(), power_neg=float(np.mean(hi < 0)),
                                         two_sided=float(np.mean((hi < 0) | (lo > 0))), mean_se_pct=float(np.mean(pc(S[v])))))
                rows_out.append(dict(site=site, model=mod, method="emp", d_true_pct=100 * d, n_windows=len(E),
                                     mean_est_pct=pc(E).mean(), power_neg=float(np.mean(E < thr)),
                                     two_sided=float(np.mean((E < thr2[0]) | (E > thr2[1]))), mean_se_pct=float("nan")))
            resid.append(dict(site=site, model=mod, null_sd_pct=float(pc(e0).std(ddof=1)), emp_thr_pct=float(pc(thr)),
                              n_windows=len(e0)))
        print("done", site, len(ws), flush=True)
    A.write_csv(os.path.join(out, "detect_calibration_injection.csv"), rows_out)
    A.write_csv(os.path.join(out, "detect_calibration_2020.csv"), true2020)
    A.write_csv(os.path.join(out, "detect_calibration_null.csv"), resid)
    A.write_csv(os.path.join(out, "detect_calibration_null_windows.csv"), nullrows)
main()
