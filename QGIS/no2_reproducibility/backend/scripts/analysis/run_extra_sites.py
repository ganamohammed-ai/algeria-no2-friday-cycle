#!/usr/bin/env python3
r"""Build the daily NO2 + ERA5 table for the extra sites of plan v2 (Analysis D: ctl_x01..04; later C) without touching
the existing outputs. Uses analyze_no2_ventilation.py unchanged, only registers the new site coordinates.
Run from the project root:  python scripts\\analysis\\run_extra_sites.py <set>
  set = D   -> tifs in data\\01_tropomi_L3_daily_geotiff\/control_extra, ERA5 in data\\03_era5\/control_extra, out outputs\/control_extra
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_no2_ventilation as A

SETS = {
    "C1": dict(sites={},  # the six study sites are already registered in analyze_no2_ventilation.py
               names="algiers,oran,arzew,constantine,setif,hassi_messaoud",
               tif=r"data/07_replication/NO2_GEE_2026", era5=r"data/07_replication/ERA5_GEE_2026",
               out=r"outputs/study_2026", years="2026"),
    "C2": dict(sites={"annaba": (36.90, 7.77), "tlemcen": (34.88, -1.32), "biskra": (34.85, 5.73), "djelfa": (34.67, 3.26),
                      "chlef": (36.17, 1.33), "tebessa": (35.40, 8.12), "tiaret": (35.37, 1.32)},
               tif=r"data/07_replication/NO2_GEE_CITIES", era5=r"data/07_replication/ERA5_GEE_CITIES",
               out=r"outputs/cities", years="2019,2020,2021,2022,2023,2024,2025,2026"),
    "D": dict(sites={"ctl_x01": (35.75, 2.25), "ctl_x02": (35.75, 3.5), "ctl_x03": (36.0, 8.0), "ctl_x04": (35.25, 7.0)},
              tif=r"data/01_tropomi_L3_daily_geotiff/control_extra", era5=r"data/03_era5/control_extra",
              out=r"outputs/control_extra", years="2019,2020,2021,2022,2023,2024,2025"),
}
import argparse
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('set',choices=['C1','C2','D']);ap.add_argument('--input-root',required=True);ap.add_argument('--out-root',required=True);args=ap.parse_args()
s=dict(SETS[args.set]);s['tif']=str(Path(args.input_root)/s['tif']);s['era5']=str(Path(args.input_root)/s['era5']);s['out']=str(Path(args.out_root)/args.set)

A.SITES.update(s["sites"])
names = s.get("names") or ",".join(s["sites"])
if "2026" in s["years"]:   # the original script stops at 31 Dec 2025; the 2026 exports end on 14 Nov 2026
    import datetime as _dt
    A.LAST_DAY = _dt.date(2026, 11, 14)
for stage in ("tables", "merge"):
    sys.argv = ["x", "--stage", stage, "--tif-dir", s["tif"], "--era5-dir", s["era5"], "--out", s["out"],
                "--sites", names, "--years", s["years"]]
    A.main()
print("done ->", s["out"])
