#!/usr/bin/env python3
"""Robustness checks for the NO2-ventilation regression (VC10): outliers, dust (UVAI), background ring, seasons.
Reads  <out>/no2_ventilation_daily.csv  (stage merge)  and  <aai-dir>/S5P_AAI2_<year>.csv.
Writes <out>/sensitivity_outliers.csv, sensitivity_dust.csv, sensitivity_season_ring.csv.
  python sensitivity_analyses.py --out outputs --aai-dir data/04_dust_aai
"""
import argparse, csv, glob, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_no2_ventilation as A

SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_rows(path):
    rows = []
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        d = {k: (v if k in ("site", "date") else fnum(v)) for k, v in r.items()}
        rows.append(d)
    return rows


def load_aai(aai_dir):
    out = {}
    for fn in sorted(glob.glob(os.path.join(aai_dir, "S5P_AAI2_*.csv"))):
        for r in csv.DictReader(open(fn, newline="", encoding="utf-8")):
            out[(r["site"], r["date"][:10])] = (fnum(r.get("uvai_inner10km")), fnum(r.get("uvai_box")))
    return out


def row(site, label, fit, key="variant"):
    if fit is None:
        return None
    return {"site": site, key: label, "n": fit["n"], "slope": round(fit["slope_lnNO2_per_lnVC"], 4), "ci95_low": round(fit["ci95_low"], 4),
            "ci95_high": round(fit["ci95_high"], 4), "partial_R2": round(fit["partial_R2"], 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs"); ap.add_argument("--aai-dir", default="data/04_dust_aai")
    ap.add_argument("--sites", nargs="+", default=["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"])
    a = ap.parse_args()
    rows = load_rows(os.path.join(a.out, "no2_ventilation_daily.csv"))
    aai = load_aai(a.aai_dir)
    for r in rows:
        r["uvai10"], r["uvaibox"] = aai.get((r["site"], r["date"]), (float("nan"), float("nan")))
        r["lnring"] = math.log(r["no2_ring"]) if r.get("no2_ring", float("nan")) > 0 else float("nan")
        r["month"] = int(r["date"][5:7])
    outl, dust, sr = [], [], []
    for s in a.sites:
        R = [r for r in rows if r["site"] == s]
        usable = [r for r in R if r["usable"] == 1 and np.isfinite(r["no2_inner"]) and np.isfinite(r["vc10"])]
        p99 = float(np.percentile([r["no2_inner"] for r in usable], 99)) if usable else float("nan")
        outl += [row(s, "all days", A.fit_site(R, "vc10"), "sample"),
                 row(s, "excl NO2>300", A.fit_site([r for r in R if r["no2_inner"] <= 300], "vc10"), "sample"),
                 row(s, "excl top 1% NO2", A.fit_site([r for r in R if r["no2_inner"] <= p99], "vc10"), "sample")]
        dust += [row(s, "base", A.fit_site(R, "vc10")),
                 row(s, "excl UVAI(10km)>1", A.fit_site([r for r in R if not r["uvai10"] > 1], "vc10")),
                 row(s, "UVAI as control", A.fit_site(R, "vc10", extra_keys=("t2m", "uvai10"))),
                 row(s, "excl UVAI(box)>1", A.fit_site([r for r in R if not r["uvaibox"] > 1], "vc10"))]
        for nm, ms in SEASONS.items():
            sr.append(row(s, "season " + nm, A.fit_site([r for r in R if r["month"] in ms], "vc10", with_month=False)))
        sr += [row(s, "base", A.fit_site(R, "vc10")), row(s, "control ln(ring)", A.fit_site(R, "vc10", extra_keys=("t2m", "lnring")))]
    for name, data in (("sensitivity_outliers.csv", outl), ("sensitivity_dust.csv", dust), ("sensitivity_season_ring.csv", sr)):
        data = [d for d in data if d]
        A.write_csv(os.path.join(a.out, name), data)
        print("wrote", name, len(data), "rows")


if __name__ == "__main__":
    main()
