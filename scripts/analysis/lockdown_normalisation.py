#!/usr/bin/env python3
"""Weather-normalised NO2 change during the 2020 containment, with placebo and control checks.

Follows analysis_plan_lockdown_normalisation.md (written before any result was computed).

Usage:
  python lockdown_normalisation.py --table no2_ventilation_daily_v2.csv \
         --controls control_no2_ventilation_daily.csv --out outputs
Writes: lockdown_effect.csv, lockdown_placebo.csv, lockdown_ml.csv, trend_version_steps.csv
"""
import argparse
import datetime as dt
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze_no2_ventilation as A  # noqa: E402

STUDY = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"]
CONTROLS = ["control_tell", "control_highland"]
D = dt.date.fromisoformat
SAMPLE_A = (D("2019-08-07"), D("2021-06-30"))
WINDOWS = {
    "primary_0323_0613": (D("2020-03-23"), D("2020-06-13")),
    "from_0329_0613": (D("2020-03-29"), D("2020-06-13")),
    "long_0323_0728": (D("2020-03-23"), D("2020-07-28")),
}
RAMADAN = [(D("2020-04-24"), D("2020-05-23")), (D("2021-04-13"), D("2021-05-12"))]
VERSION_STEPS = [D("2021-07-01"), D("2021-11-14"), D("2022-07-17"), D("2023-03-12"),
                 D("2023-11-26"), D("2024-09-08")]
PLACEBO_STEP_DAYS = 7


def prep(rows, lo=None, hi=None):
    keys = ["no2_inner", "blh", "wspd10", "t2m"]
    rr = [r for r in rows if r["usable"] == 1 and all(np.isfinite([r[k] for k in keys]))
          and r["blh"] > 0 and r["wspd10"] > 0 and r["no2_inner"] > 0]
    rr.sort(key=lambda r: r["date"])
    if lo is not None:
        rr = [r for r in rr if lo <= D(r["date"]) <= hi]
    return rr


def ramadan(dates):
    return np.array([1.0 if any(a <= D(d) <= b for a, b in RAMADAN) else 0.0 for d in dates])


def in_window(dates, win):
    return np.array([1.0 if win[0] <= D(d) <= win[1] else 0.0 for d in dates])


def fit(rr, win, adjust, placebo_win=None):
    """OLS with Newey-West errors. The effect of interest is always design column 1.
    Returns (L, se, n_days_in_target_window)."""
    dates = [r["date"] for r in rr]
    y = np.log([r["no2_inner"] for r in rr])
    lb = np.log([r["blh"] for r in rr])
    lw = np.log([r["wspd10"] for r in rr])
    t2 = np.array([r["t2m"] for r in rr])
    real = in_window(dates, win)
    target = in_window(dates, placebo_win) if placebo_win is not None else real
    extra = []
    if adjust:
        extra += [lw, lb, t2]
    extra.append(ramadan(dates))
    if placebo_win is not None:
        extra.append(real)  # the true lockdown stays in the model
    X, _ = A.design(dates, target, extra, with_month=True)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    V = A.newey_west(X, res, A.NW_LAG)
    return float(beta[1]), float(math.sqrt(V[1, 1])), int(target.sum())


def pct(l):
    return 100.0 * (math.exp(l) - 1.0)


def lockdown_table(sitedata, label):
    out = []
    for site, rows in sitedata.items():
        rr = prep(rows, *SAMPLE_A)
        for wname, win in WINDOWS.items():
            n_in = sum(1 for r in rr if win[0] <= D(r["date"]) <= win[1])
            if n_in < 15:
                out.append(dict(group=label, site=site, window=wname, n_days_window=n_in))
                continue
            lr, sr, _ = fit(rr, win, False)
            la, sa, _ = fit(rr, win, True)
            out.append(dict(group=label, site=site, window=wname, n_days_window=n_in, n_total=len(rr),
                            raw_pct=pct(lr), raw_lo=pct(lr - 1.96 * sr), raw_hi=pct(lr + 1.96 * sr),
                            adj_pct=pct(la), adj_lo=pct(la - 1.96 * sa), adj_hi=pct(la + 1.96 * sa),
                            weather_part_pct_points=pct(lr) - pct(la)))
    return out


def placebo_table(sitedata, label):
    out = []
    win = WINDOWS["primary_0323_0613"]
    length = (win[1] - win[0]).days
    for site, rows in sitedata.items():
        rr = prep(rows, *SAMPLE_A)
        raw, adj = [], []
        start = SAMPLE_A[0]
        while start + dt.timedelta(days=length) <= SAMPLE_A[1]:
            end = start + dt.timedelta(days=length)
            overlap = not (end < win[0] or start > win[1])
            n_in = sum(1 for r in rr if start <= D(r["date"]) <= end)
            if not overlap and n_in >= 15:
                try:
                    raw.append(fit(rr, win, False, placebo_win=(start, end))[0])
                    adj.append(fit(rr, win, True, placebo_win=(start, end))[0])
                except np.linalg.LinAlgError:
                    pass
            start += dt.timedelta(days=PLACEBO_STEP_DAYS)
        if len(raw) < 10:
            out.append(dict(group=label, site=site, n_placebo=len(raw)))
            continue
        lr, la = fit(rr, win, False)[0], fit(rr, win, True)[0]
        raw, adj = np.array(raw), np.array(adj)
        out.append(dict(group=label, site=site, n_placebo=len(raw),
                        sd_raw_pct=float(np.std(100 * (np.exp(raw) - 1), ddof=1)),
                        sd_adj_pct=float(np.std(100 * (np.exp(adj) - 1), ddof=1)),
                        rank_raw=float(np.mean(np.abs(raw) >= abs(lr))),
                        rank_adj=float(np.mean(np.abs(adj) >= abs(la)))))
    return out


def ml_table(sitedata, label, seed=0):
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.metrics import r2_score
    out = []
    win = WINDOWS["primary_0323_0613"]
    rng = np.random.default_rng(seed)
    for site, rows in sitedata.items():
        rr = prep(rows, *SAMPLE_A)
        dates = [D(r["date"]) for r in rr]
        y = np.log([r["no2_inner"] for r in rr])
        feats = np.column_stack([
            np.log([r["wspd10"] for r in rr]), np.log([r["blh"] for r in rr]),
            [r["t2m"] for r in rr], [r["u10"] for r in rr], [r["v10"] for r in rr],
            [d.timetuple().tm_yday for d in dates], [d.weekday() for d in dates],
            [(d - D("2019-01-01")).days / 365.25 for d in dates]])
        inwin = np.array([win[0] <= d <= win[1] for d in dates])
        if inwin.sum() < 15:
            out.append(dict(group=label, site=site, n_days_window=int(inwin.sum())))
            continue
        train = ~inwin
        mk = lambda: HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_depth=4,
                                                    min_samples_leaf=20, random_state=0)
        m = mk().fit(feats[train], y[train])
        resid_win = y[inwin] - m.predict(feats[inwin])
        # time-blocked CV (6 blocks) on training days
        idx = np.where(train)[0]
        blocks = np.array_split(idx, 6)
        pred = np.full(len(y), np.nan)
        for b in blocks:
            tr = np.setdiff1d(idx, b)
            pred[b] = mk().fit(feats[tr], y[tr]).predict(feats[b])
        r2 = r2_score(y[idx], pred[idx])
        resid_train = (y[idx] - pred[idx])
        L = float(resid_win.mean())
        # block bootstrap of out-of-sample residuals (7-day blocks) for the mean over n_win days
        nwin = int(inwin.sum())
        means = []
        nb = max(1, nwin // 7)
        starts = np.arange(0, len(resid_train) - 7)
        for _ in range(500):
            s = rng.choice(starts, size=nb)
            means.append(np.concatenate([resid_train[i:i + 7] for i in s]).mean())
        lo, hi = np.percentile(means, [2.5, 97.5])
        out.append(dict(group=label, site=site, n_days_window=nwin, cv_r2=float(r2),
                        ml_pct=pct(L), null_lo_pct=pct(lo), null_hi_pct=pct(hi)))
    return out


def trend_table(sitedata, label):
    """Full period 2019-08-07..2025-12-31; trend %/yr raw vs adjusted, with/without version steps."""
    out = []
    for site, rows in sitedata.items():
        rr = prep(rows)
        dates = [r["date"] for r in rr]
        y = np.log([r["no2_inner"] for r in rr])
        lb = np.log([r["blh"] for r in rr])
        lw = np.log([r["wspd10"] for r in rr])
        t2 = np.array([r["t2m"] for r in rr])
        steps = [np.array([1.0 if D(d) >= s else 0.0 for d in dates]) for s in VERSION_STEPS]
        covid = in_window(dates, WINDOWS["long_0323_0728"])
        rm = ramadan(dates)
        res_row = dict(group=label, site=site, n=len(rr))
        for adjust in (False, True):
            for with_steps in (False, True):
                extra = []
                if adjust:
                    extra += [lw, lb, t2]
                extra.append(rm)
                if with_steps:
                    extra += steps
                X, names = A.design(dates, covid, extra, with_month=True)
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                res = y - X @ beta
                V = A.newey_west(X, res, A.NW_LAG)
                ti = names.index("trend_per_year")
                b_, s_ = beta[ti], math.sqrt(V[ti, ti])
                tag = ("adj" if adjust else "raw") + ("_steps" if with_steps else "_nosteps")
                res_row[tag + "_pct_per_yr"] = pct(b_)
                res_row[tag + "_lo"] = pct(b_ - 1.96 * s_)
                res_row[tag + "_hi"] = pct(b_ + 1.96 * s_)
        out.append(res_row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--out", default=".")
    ap.add_argument("--skip-ml", action="store_true", help="Replay OLS, placebo and trend tables without the expensive ML fits")
    a = ap.parse_args()
    study_rows = A.read_csv(a.table)
    ctrl_rows = A.read_csv(a.controls)
    groups = {
        "study": {s: [r for r in study_rows if r["site"] == s] for s in STUDY},
        "control": {s: [r for r in ctrl_rows if r["site"] == s] for s in CONTROLS},
    }
    L, P, M, T = [], [], [], []
    for label, sd in groups.items():
        sd = {k: v for k, v in sd.items() if v}
        L += lockdown_table(sd, label)
        P += placebo_table(sd, label)
        if not a.skip_ml:
            M += ml_table(sd, label)
        T += trend_table(sd, label)
    for name, rows in (("lockdown_effect.csv", L), ("lockdown_placebo.csv", P),
                       ("lockdown_ml.csv", M), ("trend_version_steps.csv", T)):
        if name == "lockdown_ml.csv" and a.skip_ml:
            continue
        A.write_csv(os.path.join(a.out, name), rows)
        print("wrote", name, len(rows))


if __name__ == "__main__":
    main()
