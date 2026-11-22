#!/usr/bin/env python3
"""Which TROPOMI NO2 processor versions are in the Level-2 clips? (reads progress.csv only; no download)

progress.csv is written by scripts\\download\\tropomi_download.py. Each row is one orbit file with its
class (RPRO = reprocessed, OFFL = operational) and its 6-digit processor version (020400 = v2.4.0).
Output: a table by class and version (first/last day, number of files and of days) and a flag
for the versions that KNMI/TEMIS treat as equivalent to the v2.4.0 reprocessing (v2.4.0 and v2.5.0).

Usage (from the project folder):
  python scripts\\analysis\\version_audit.py --progress data\\02_tropomi_L2_orbit_clips\\progress.csv --out outputs
"""
import argparse, csv, collections, os

# v2.5.0 only fixes the OFFL qa_value over snow/ice and is described by KNMI/TEMIS as consistent with RPRO v2.4.0.
SAME_AS_V240 = {"020400", "020500"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--progress", required=True); ap.add_argument("--out", default=".")
    a = ap.parse_args()
    g = collections.defaultdict(lambda: {"files": 0, "days": set()})
    n_bad = 0
    with open(a.progress, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["status"] != "done":
                n_bad += 1; continue
            key = (r["mode"], str(r["processor"]).zfill(6))
            g[key]["files"] += 1; g[key]["days"].add(r["date"])
    rows = []
    for (mode, proc), v in sorted(g.items(), key=lambda kv: min(kv[1]["days"])):
        d = sorted(v["days"])
        rows.append({"class": mode, "processor": "v%d.%d.%d" % (int(proc[:2]), int(proc[2:4]), int(proc[4:])),
                     "first_day": d[0], "last_day": d[-1], "files": v["files"], "days": len(d),
                     "same_as_v2.4.0": "yes" if proc in SAME_AS_V240 else "no"})
    print("%-5s %-9s %-9s %-9s %6s %6s  %s" % ("class", "version", "first", "last", "files", "days", "same as v2.4.0"))
    for r in rows:
        print("%-5s %-9s %-9s %-9s %6d %6d  %s" % (r["class"], r["processor"], r["first_day"], r["last_day"],
                                                  r["files"], r["days"], r["same_as_v2.4.0"]))
    same = [r for r in rows if r["same_as_v2.4.0"] == "yes"]
    ndays = sum(r["days"] for r in same); last = max(r["last_day"] for r in same)
    print("\nFiles not 'done' (skipped): %d" % n_bad)
    print("Days in v2.4.0-equivalent versions: %d, up to %s" % (ndays, last))
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, "l2_version_audit.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("wrote", out)

if __name__ == "__main__":
    main()
