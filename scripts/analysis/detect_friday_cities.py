#!/usr/bin/env python3
"""Analysis C, holdout 2 of analysis_plan_v2_detectability.md: the Friday model of friday_cycle.py, unchanged, on seven new cities,
7 Aug 2019 - 14 Nov 2026 (the cities were never analysed before). Only change: Ramadan 2026 (approx. 18 Feb - 19 Mar) added.
Prediction (written 20 Nov 2026): negative effect in all cities; Annaba (heavy industry) smaller in magnitude than the six traffic/residential cities."""
import sys, csv
import friday_cycle as F
A = F.A
F.RAM.append((F.D("2026-02-18"), F.D("2026-03-19")))
import argparse,os
ap=argparse.ArgumentParser(); ap.add_argument("table"); ap.add_argument("--out",default="."); args=ap.parse_args(); os.makedirs(args.out,exist_ok=True)
rows = A.read_csv(args.table)
CITIES = ["annaba", "tlemcen", "biskra", "djelfa", "chlef", "tebessa", "tiaret"]
out = []
print("city        n   Fri ring-controlled (95% CI)      weather-adjusted        Sat ring")
for s in CITIES:
    sub = [r for r in rows if r["site"] == s]
    res = {}
    for adj, tag in ((1, "adj"), (2, "ring")):
        n, o = F.fit(sub, adj, "fs"); res[tag] = o["Fri"]; res[tag+"_sat"] = o["Sat"]; res["n"] = n
    print("%-10s %4d  %6.1f (%6.1f..%6.1f)   %6.1f (%6.1f..%6.1f)   %6.1f (%6.1f..%6.1f)" % (s, res["n"], *res["ring"], *res["adj"], *res["ring_sat"]))
    out.append(dict(site=s, n=res["n"], fri_ring_pct=res["ring"][0], fri_ring_lo=res["ring"][1], fri_ring_hi=res["ring"][2],
                    fri_adj_pct=res["adj"][0], fri_adj_lo=res["adj"][1], fri_adj_hi=res["adj"][2],
                    sat_ring_pct=res["ring_sat"][0], sat_ring_lo=res["ring_sat"][1], sat_ring_hi=res["ring_sat"][2]))
A.write_csv(os.path.join(args.out,"detect_friday_cities.csv"), out)
