#!/usr/bin/env python3
"""Analysis A of analysis_plan_v2_detectability.md: synthetic injection into placebo windows.

Reuses prep/fit/WINDOWS/SAMPLE_A of lockdown_normalisation.py unchanged. The injection multiplies
no2_inner on the days of a placebo window by (1+d).

Usage: python detect_injection.py --table T --controls C --out OUT [--gb]
"""
import argparse, datetime as dt, os, sys
import numpy as np
import lockdown_normalisation as L  # noqa
import analyze_no2_ventilation as A  # noqa

D = L.D
DELTAS = [0.0, -0.05, -0.10, -0.15, -0.20, -0.30]
WIN = L.WINDOWS["primary_0323_0613"]
LENGTH = (WIN[1] - WIN[0]).days


def windows(rr):
    out, start = [], L.SAMPLE_A[0]
    while start + dt.timedelta(days=LENGTH) <= L.SAMPLE_A[1]:
        end = start + dt.timedelta(days=LENGTH)
        overlap = not (end < WIN[0] or start > WIN[1])
        n_in = sum(1 for r in rr if start <= D(r["date"]) <= end)
        if not overlap and n_in >= 15:
            out.append((start, end))
        start += dt.timedelta(days=L.PLACEBO_STEP_DAYS)
    return out


def inject(rr, w, d):
    new = []
    for r in rr:
        if w[0] <= D(r["date"]) <= w[1]:
            r = dict(r); r["no2_inner"] = r["no2_inner"] * (1.0 + d)
        new.append(r)
    return new


def gb_shift(rr, w):
    """Mean observed-minus-predicted ln NO2 over the window (v1 features), no CV, no bootstrap."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    dates = [D(r["date"]) for r in rr]
    y = np.log([r["no2_inner"] for r in rr])
    f = np.column_stack([
        np.log([r["wspd10"] for r in rr]), np.log([r["blh"] for r in rr]), [r["t2m"] for r in rr],
        [r["u10"] for r in rr], [r["v10"] for r in rr], [d.timetuple().tm_yday for d in dates],
        [d.weekday() for d in dates], [(d - D("2019-01-01")).days / 365.25 for d in dates]])
    inw = np.array([w[0] <= d <= w[1] for d in dates])
    real = np.array([WIN[0] <= d <= WIN[1] for d in dates])
    train = ~inw & ~real
    m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_depth=4,
                                      min_samples_leaf=20, random_state=0).fit(f[train], y[train])
    return float((y[inw] - m.predict(f[inw])).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table"); ap.add_argument("--controls"); ap.add_argument("--out", default=".")
    ap.add_argument("--gb", action="store_true"); ap.add_argument("--sites", default="")
    ap.add_argument("--check", action="store_true", help="reproduce v1 placebo SDs only")
    a = ap.parse_args()
    st, ct = A.read_csv(a.table), A.read_csv(a.controls)
    groups = [(s, "study", [r for r in st if r["site"] == s]) for s in L.STUDY] + \
             [(s, "control", [r for r in ct if r["site"] == s]) for s in L.CONTROLS]
    if a.sites:
        groups = [g for g in groups if g[0] in a.sites.split(",")]
    res, gbres = [], []
    for site, grp, rows in groups:
        rr = L.prep(rows, *L.SAMPLE_A)
        ws = windows(rr)
        deltas = [0.0] if a.check else DELTAS
        for d in deltas:
            for adjust in (False, True):
                est, lo, hi = [], [], []
                for w in ws:
                    rr2 = inject(rr, w, d)
                    try:
                        l, s, _ = L.fit(rr2, WIN, adjust, placebo_win=w)
                    except np.linalg.LinAlgError:
                        continue
                    est.append(l); lo.append(l - 1.96 * s); hi.append(l + 1.96 * s)
                est, lo, hi = map(np.array, (est, lo, hi))
                pc = lambda x: 100 * (np.exp(x) - 1)
                res.append(dict(site=site, group=grp, model="adjusted" if adjust else "unadjusted",
                                d_true_pct=100 * d, n_windows=len(est), mean_est_pct=float(pc(est).mean()),
                                sd_est_pct=float(pc(est).std(ddof=1)),
                                bias_pct_points=float(pc(est).mean() - 100 * d),
                                power=float(np.mean((hi < 0) | (lo > 0))),
                                power_neg=float(np.mean(hi < 0))))
        if a.gb and not a.check:
            for d in DELTAS:
                sh = np.array([gb_shift(inject(rr, w, d), w) for w in ws])
                gbres.append(dict(site=site, group=grp, d_true_pct=100 * d, n_windows=len(sh),
                                  mean_est_pct=float(100 * (np.exp(sh) - 1).mean()), _raw=sh))
        print("done", site, len(ws), flush=True)
    A.write_csv(os.path.join(a.out, "detect_injection_ols.csv"), res)
    if gbres:
        # empirical null: 2.5/97.5 percentile of GB shifts at d=0 in the same site
        out = []
        for g in gbres:
            if g["d_true_pct"] == 0.0:
                null = g["_raw"]
        for site in {g["site"] for g in gbres}:
            base = [g for g in gbres if g["site"] == site and g["d_true_pct"] == 0.0][0]["_raw"]
            lo_, hi_ = np.percentile(base, [2.5, 97.5])
            for g in [g for g in gbres if g["site"] == site]:
                out.append(dict(site=site, group=g["group"], d_true_pct=g["d_true_pct"],
                                n_windows=g["n_windows"], mean_est_pct=g["mean_est_pct"],
                                power_empirical=float(np.mean(g["_raw"] < lo_)),
                                false_pos_rate_neg=float(np.mean(base < lo_))))
        A.write_csv(os.path.join(a.out, "detect_injection_gb.csv"), out)


if __name__ == "__main__":
    main()
