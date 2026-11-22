#!/usr/bin/env python3
"""Analysis C, holdout 1 of analysis_plan_v2_detectability.md: the Friday model of friday_cycle.py, unchanged, on 2026 data only.
Only change: Ramadan 2026 (approx. 18 Feb - 19 Mar 2026) is added to the Ramadan dates. Prediction (written 20 Nov 2026, before the data were opened):
negative Friday effect with 95% interval excluding zero at Algiers, Oran, Constantine, Setif; Arzew and Hassi Messaoud smaller in magnitude than those four."""
import sys, csv, math
import friday_cycle as F
A = F.A
F.RAM.append((F.D("2026-02-18"), F.D("2026-03-19")))
import argparse,os
ap=argparse.ArgumentParser(); ap.add_argument("table"); ap.add_argument("reference"); ap.add_argument("--out",default="."); args=ap.parse_args(); os.makedirs(args.out,exist_ok=True)
rows = A.read_csv(args.table)
old = {r["site"]: r for r in csv.DictReader(open(args.reference))}
out = []
print("site            n   Fri ring-controlled 2026 (95% CI)      adjusted-only          | 2019-2025 ring")
for s in F.STUDY:
    sub = [r for r in rows if r["site"] == s]
    res = {}
    for adj, tag in ((1, "adj"), (2, "ring")):
        n, o = F.fit(sub, adj, "fs"); res[tag] = o["Fri"]; res["n"] = n
    o0 = old[s]
    print("%-15s %4d  %6.1f (%6.1f..%6.1f)   %6.1f (%6.1f..%6.1f) | %6.1f" % (s, res["n"], *res["ring"], *res["adj"], float(o0["Fri_ring_pct"])))
    out.append(dict(site=s, n=res["n"], fri_ring_pct=res["ring"][0], fri_ring_lo=res["ring"][1], fri_ring_hi=res["ring"][2],
                    fri_adj_pct=res["adj"][0], fri_adj_lo=res["adj"][1], fri_adj_hi=res["adj"][2], fri_ring_2019_2025=float(o0["Fri_ring_pct"])))
A.write_csv(os.path.join(args.out,"detect_friday_holdout_2026.csv"), out)
