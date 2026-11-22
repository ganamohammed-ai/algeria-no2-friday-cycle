#!/usr/bin/env python3
"""Daily NO2 table built from the Level-2 orbit clips, restricted to ONE processor version.

Why: the Google Earth Engine Level-3 collection mixes processor versions (v1.3 -> v1.4.0 on 2020-11-29,
v2.2.0 on 2021-07-01, ...). The clips downloaded from Copernicus hold the reprocessed RPRO v2.4.0 record
(2019-08-07 .. 2022-07-25) followed by OFFL v2.4.0 and v2.5.0 (to 2023-11-25), which KNMI/TEMIS treat as
the same retrieval. Using only those files gives one version for 1566 days.

For every site and day it computes, from pixels with qa_value >= --qa:
  no2_inner  mean tropospheric column in the 10 km disc (umol m-2, like the Level-3 table)
  no2_ring   mean in the 50-100 km ring
A day is 'usable' if, in BOTH areas, valid pixels are at least 30% of the stored pixels with a retrieval
(same 30% rule as the Level-3 table; the denominator differs because the clips keep only pixels with a retrieval)
and there are at least 3 valid pixels in the disc.
The ERA5 columns are copied from the existing daily table (--table), so weather is identical to the main analysis.

With the controls: --clips CLIPS_STUDY CLIPS_CONTROLS --table STUDY_TABLE CONTROL_TABLE.

Usage (from the project folder):
  python scripts\\analysis\\l2_daily_table.py --clips data\\02_tropomi_L2_orbit_clips --table outputs\\no2_ventilation_daily_v2.csv --qa 0.75 --out outputs --tag q75
"""
import argparse, collections, csv, glob, math, os
import numpy as np

SITES = {"algiers": (36.75, 3.06), "oran": (35.70, -0.63), "arzew": (35.85, -0.31),
         "constantine": (36.37, 6.61), "setif": (36.19, 5.41), "hassi_messaoud": (31.68, 6.07),
         "control_tell": (35.0, 5.25), "control_highland": (34.5, -1.0)}
INNER_KM, RING_KM = 10.0, (50.0, 100.0)
MIN_FRAC, MIN_INNER_PIX = 0.30, 3
KEEP_VERSIONS = {"020400", "020500"}
R_EARTH = 6371.0


def one_file(z, lat0, lon0, qa_min):
    lat = z["latitude"].astype(float); lon = z["longitude"].astype(float)
    no2 = z["no2_trop_mol_m2"].astype(float) * 1e6
    qa = z["qa_value"].astype(float)
    dx = np.radians(lon - lon0) * R_EARTH * math.cos(math.radians(lat0))
    dy = np.radians(lat - lat0) * R_EARTH
    d = np.hypot(dx, dy)
    have = np.isfinite(no2)
    val = have & (qa >= qa_min)
    inner, ring = d <= INNER_KM, (d >= RING_KM[0]) & (d <= RING_KM[1])
    return {"in_sum": no2[inner & val].sum(), "in_n": int((inner & val).sum()), "in_tot": int((inner & have).sum()),
            "rg_sum": no2[ring & val].sum(), "rg_n": int((ring & val).sum()), "rg_tot": int((ring & have).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clips", required=True, nargs="+", help="one or more clip folders (study sites, controls)")
    ap.add_argument("--table", required=True, nargs="+", help="one or more daily tables with the ERA5 columns")
    ap.add_argument("--out", default="."); ap.add_argument("--qa", type=float, default=0.75)
    ap.add_argument("--tag", default="q75")
    a = ap.parse_args()
    base, cols = {}, []
    for tpath in a.table:
        with open(tpath, newline="", encoding="utf-8") as f:
            rd = csv.DictReader(f)
            cols += [c for c in rd.fieldnames if c not in cols]
            for r in rd:
                base[(r["site"], r["date"])] = r
    out_rows, skipped = [], collections.Counter()
    for site, (la, lo) in SITES.items():
        files = sorted(f for c in a.clips for f in glob.glob(os.path.join(c, site, "*.npz")))
        agg = {}
        for i, fp in enumerate(files):
            d8 = os.path.basename(fp)[:8]; date = "%s-%s-%s" % (d8[:4], d8[4:6], d8[6:8])
            try:
                z = np.load(fp, allow_pickle=True)
                ver = str(z["processor_version"]) if "processor_version" in z.files else ""
                if ver not in KEEP_VERSIONS:
                    skipped["other_version"] += 1; continue
                res = one_file(z, la, lo, a.qa)
            except Exception as e:
                skipped["error"] += 1; print("skip", fp, e); continue
            c = agg.setdefault(date, collections.Counter())
            for k, v in res.items():
                c[k] += v
            if (i + 1) % 300 == 0:
                print(site, i + 1, "/", len(files), flush=True)
        n_days = 0
        for date in sorted(agg):
            c = agg[date]; b = base.get((site, date))
            if b is None:
                skipped["no_era5_row"] += 1; continue
            if c["in_n"] == 0 or c["rg_n"] == 0:
                skipped["no_valid_pixels"] += 1; continue
            ok = (c["in_n"] >= MIN_FRAC * c["in_tot"] and c["rg_n"] >= MIN_FRAC * c["rg_tot"] and c["in_n"] >= MIN_INNER_PIX)
            ni, nr = c["in_sum"] / c["in_n"], c["rg_sum"] / c["rg_n"]
            r = dict(b)
            r.update({"no2_inner": ni, "inner_n": c["in_n"], "inner_total": c["in_tot"], "no2_ring": nr,
                      "ring_n": c["rg_n"], "ring_total": c["rg_tot"], "no2_box": "", "box_n": "",
                      "no2_delta": ni - nr, "usable": 1 if ok else 0})
            out_rows.append(r); n_days += 1
        print(site, "days with L2 v2.4.0-equivalent data:", n_days)
    path = os.path.join(a.out, "l2_v240_daily_%s.csv" % a.tag)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, restval=""); w.writeheader(); w.writerows(out_rows)
    print("wrote", path, "rows", len(out_rows), "usable", sum(int(r["usable"]) for r in out_rows), "skipped", dict(skipped))


if __name__ == "__main__":
    main()
