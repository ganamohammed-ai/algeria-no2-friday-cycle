#!/usr/bin/env python3
"""Fixed rule for choosing the control (background) sites. Uses ONLY static descriptors, never NO2.
Input : Control_candidates.csv (from Control_site_screen.js)
Output: control_sites_selected.csv  + a printed audit table.

Rule (written before the candidate table was seen):
 1. Eligible = built-up share within 10 km <= 0.5%  AND  within 50 km <= 2%
               AND population within 10 km <= 5,000
               AND night-light within 10 km <= 1.0 nW cm-2 sr-1
               AND water share within 10 km <= 5%
               AND distance to the nearest of the six study sites >= 100 km.
 2. Two controls, one per elevation stratum: 'tell' (mean elevation 0-600 m) and 'highland' (>= 800 m).
    In each stratum take the eligible candidate with the lowest night-light within 10 km;
    ties (within 0.05) are broken by the lowest built-up share, then by the smaller elevation s.d.
 3. If a stratum has no eligible candidate, relax in this order, one step at a time, until one exists:
    distance >= 80 km, night-light <= 2.0, population <= 10,000, built-up 50 km <= 5%.

DEVIATION NOTE (made after the static candidate table was seen, before any NO2 or wind data for the candidates was opened):
 The literal rule above selected points at 33.5 N, about 270 km from the nearest study site, in the Sahara depressions
 (one outside Algeria), which share neither climate nor synoptic regime with the northern sites. Run with --amended, the rule
 adds two NO2-blind conditions: latitude >= 34.5 N and distance to the nearest study site <= 150 km (the 100 km minimum stays).
 Both outputs are kept; the amended one is the control used in the paper, and the deviation is reported.
"""
import csv, sys

LADDER = [dict(), dict(dist=80), dict(dist=80, nl=2.0), dict(dist=80, nl=2.0, pop=10000), dict(dist=80, nl=2.0, pop=10000, b50=0.05)]
BASE = dict(b10=0.005, b50=0.02, pop=5000, nl=1.0, water=0.05, dist=100)


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


AMENDED = False          # set by --amended (see the deviation note in the docstring)


def pick(rows, lo, hi):
    for step in LADDER:
        p = dict(BASE); p.update(step)
        el = [r for r in rows if (not AMENDED or (r["lat"] >= 34.5 and r["dist"] <= 150)) and lo <= r["elev"] < hi and r["b10"] <= p["b10"] and r["b50"] <= p["b50"] and r["pop"] <= p["pop"]
              and r["nl"] <= p["nl"] and r["water"] <= p["water"] and r["dist"] >= p["dist"]]
        if el:
            m = min(r["nl"] for r in el)
            c = [r for r in el if r["nl"] <= m + 0.05]
            c.sort(key=lambda r: (r["b10"], r["elev_sd"]))
            return c[0], step, len(el)
    return None, None, 0


def main(path, out="control_sites_selected.csv"):
    global AMENDED
    AMENDED = "--amended" in sys.argv
    if AMENDED:
        out = "control_sites_selected_amended.csv"
    rows = []
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        d = dict(lat=f(r["lat"]), lon=f(r["lon"]), b10=f(r["built_10km"]), b50=f(r["built_50km"]), water=f(r.get("water_10km", 0) or 0),
                 pop=f(r["pop_10km"]), nl=f(r["nightlight_10km"]), elev=f(r["elev_mean_10km"]), elev_sd=f(r["elev_sd_10km"]),
                 dist=f(r["dist_nearest_site_km"]))
        if any(v != v for k, v in d.items() if k != "water"):
            continue
        if d["water"] != d["water"]:
            d["water"] = 0.0
        rows.append(d)
    res = []
    for name, lo, hi in (("tell", 0, 600), ("highland", 800, 5000)):
        r, step, n = pick(rows, lo, hi)
        print(name, "eligible:", n, "relaxation:", step, "->", r)
        if r:
            res.append(dict(stratum=name, relaxation=str(step), n_eligible=n, **r))
    if res:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
        print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Control_candidates.csv")
