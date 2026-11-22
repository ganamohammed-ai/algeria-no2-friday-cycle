#!/usr/bin/env python3
"""Download hourly ERA5 meteorology for the study sites from the Copernicus
Climate Data Store (CDS), only at the TROPOMI overpass hours (12:00 and 13:00 UTC).

One-time setup (see README_AR.txt):
  1. free account at https://cds.climate.copernicus.eu
  2. copy your Personal Access Token from your profile page into a file called
     .cdsapirc in your home folder, two lines:
         url: https://cds.climate.copernicus.eu/api
         key: <your-token>
  3. on the dataset pages "ERA5 hourly data on single levels" and
     "ERA5 hourly data on pressure levels", tick "Accept terms" once.

Examples:
    python era5_download.py --years 2024 --months 1          # small test
    python era5_download.py --years 2019-2025                # the TROPOMI period
    python era5_download.py --years 1991-2020 --per-site     # climate reference (slow, run last)

Already downloaded months are skipped, so it is safe to run again.
"""
import argparse
import calendar
import re
import sys
import time
from pathlib import Path

from sites import SITES, parse_site_list, site_box, HALF_BOX_DEG

SINGLE_VARS = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "boundary_layer_height",
    "2m_temperature",
    "surface_pressure",
    "surface_solar_radiation_downwards",
    "total_precipitation",
]
PL_VARS = ["u_component_of_wind", "v_component_of_wind", "geopotential"]
LEVELS = ["1000", "975", "950", "925", "900", "850", "500"]
HOURS = ["12:00", "13:00"]  # TROPOMI passes over Algeria at ~13:30 local = ~12:30 UTC


def parse_years(text):
    if "-" in text:
        a, b = text.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(t) for t in text.split(",")]


def area_for(site_names, half, margin=0.25):
    """[North, West, South, East] covering all given sites (CDS order)."""
    boxes = [site_box(s, half) for s in site_names]
    w = min(b[0] for b in boxes) - margin
    s = min(b[1] for b in boxes) - margin
    e = max(b[2] for b in boxes) + margin
    n = max(b[3] for b in boxes) + margin
    r = lambda v: round(v * 4) / 4  # snap to the 0.25 degree ERA5 grid
    return [r(n), r(w), r(s), r(e)]


def covered_months(out, kind, g, y):
    """Months of year y already present for this kind, from any file naming (yearly, range, monthly)."""
    got = set()
    for p in out.glob("era5_%s_%s_%d*.nc" % (kind, g, y)):
        if p.stat().st_size == 0:
            continue
        n = p.stem
        m = re.search(r"_%d_m(\d\d)-(\d\d)$" % y, n)
        if m:
            got.update(range(int(m.group(1)), int(m.group(2)) + 1))
            continue
        m = re.search(r"_%d_(\d\d)$" % y, n)
        if m:
            got.add(int(m.group(1)))
            continue
        if re.search(r"_%d$" % y, n):
            got.update(range(1, 13))
    return got


def target_name(out, kind, g, y, a, b):
    if a == 1 and b == 12:
        return out / ("era5_%s_%s_%d.nc" % (kind, g, y))
    if a == b:
        return out / ("era5_%s_%s_%d_%02d.nc" % (kind, g, y, a))
    return out / ("era5_%s_%s_%d_m%02d-%02d.nc" % (kind, g, y, a, b))


def fetch(client, out, kind, g, area, y, a, b, failed, depth=0):
    """Download months a..b of year y. Waits when the CDS queue is full, splits the request in two
    when the CDS says it is too large."""
    have = covered_months(out, kind, g, y)
    missing = [m for m in range(a, b + 1) if m not in have]
    if not missing:
        return
    a, b = missing[0], missing[-1]      # ask only for the part that is not on disk yet
    target = target_name(out, kind, g, y, a, b)
    days = ["%02d" % d for d in range(1, 32)] if b > a else ["%02d" % d for d in range(1, calendar.monthrange(y, a)[1] + 1)]
    req = {
        "product_type": ["reanalysis"],
        "variable": SINGLE_VARS if kind == "sl" else PL_VARS,
        "year": [str(y)],
        "month": ["%02d" % m for m in range(a, b + 1)],
        "day": days,
        "time": HOURS,
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    if kind == "pl":
        req["pressure_level"] = LEVELS
    name = "reanalysis-era5-single-levels" if kind == "sl" else "reanalysis-era5-pressure-levels"
    label = "%d months %02d-%02d" % (y, a, b)
    waits = 0
    while True:
        print("%s %s %s" % (time.strftime("%H:%M:%S"), name.replace("reanalysis-era5-", ""), label), flush=True)
        try:
            client.retrieve(name, req, str(target))
            return
        except Exception as e:
            msg = str(e)
            if target.exists() and target.stat().st_size == 0:
                target.unlink()
            low = msg.lower()
            if "queued requests" in low or "maximum number" in low:
                waits += 1
                if waits > 300:
                    print("    gave up waiting for a free queue slot", flush=True)
                    failed.append(target.name)
                    return
                print("    CDS queue is full for this dataset; waiting 5 minutes then trying again (%d)" % waits, flush=True)
                time.sleep(300)
                continue
            if ("cost limits" in low or "too large" in low) and b > a:
                mid = (a + b) // 2
                print("    request too large: splitting %02d-%02d into %02d-%02d and %02d-%02d" % (a, b, a, mid, mid + 1, b), flush=True)
                fetch(client, out, kind, g, area, y, a, mid, failed, depth + 1)
                fetch(client, out, kind, g, area, y, mid + 1, b, failed, depth + 1)
                return
            print("    FAILED: %s" % msg.splitlines()[0] if msg else "    FAILED", flush=True)
            failed.append(target.name)
            return


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", default="2019-2025", help="e.g. 2019-2025 or 2024 or 2020,2021")
    ap.add_argument("--months", default="1-12", help="used with --chunk month, e.g. 1-12 or 1 or 12,1,2")
    ap.add_argument("--sites", default="default")
    ap.add_argument("--half-box", type=float, default=HALF_BOX_DEG)
    ap.add_argument("--per-site", action="store_true",
                    help="one small request per site instead of one box covering all sites")
    ap.add_argument("--outdir", default="era5")
    ap.add_argument("--skip-pressure", action="store_true", help="single-level variables only")
    ap.add_argument("--only", choices=["sl", "pl"], help="download only the single-level (sl) or pressure-level (pl) files")
    ap.add_argument("--chunk", choices=["month", "quarter", "year"], default="month",
                    help="months per request. year / quarter = far fewer requests (each waits in the CDS queue "
                         "only once); a request the CDS finds too large is split automatically")
    args = ap.parse_args()

    try:
        import cdsapi
    except ImportError:
        raise SystemExit("Missing package. Run:  pip install cdsapi")

    years = parse_years(args.years)
    if "-" in args.months:
        a, b = args.months.split("-")
        months = list(range(int(a), int(b) + 1))
    else:
        months = [int(t) for t in args.months.split(",")]
    sites = parse_site_list(args.sites)
    groups = [(s, area_for([s], args.half_box)) for s in sites] if args.per_site \
        else [("all", area_for(sites, args.half_box))]

    kinds = ["sl"] if args.skip_pressure else ["sl", "pl"]
    if args.only:
        kinds = [args.only]
    if args.chunk == "year":
        spans = [(1, 12)]
    elif args.chunk == "quarter":
        spans = [(1, 3), (4, 6), (7, 9), (10, 12)]
    else:
        spans = [(m, m) for m in months]

    client = cdsapi.Client()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(g, area, y, a, b, kind) for g, area in groups for y in years for (a, b) in spans for kind in kinds]
    print("Requests to make: %d  (each waits in the CDS queue; this can take hours)" % len(jobs), flush=True)

    failed = []
    for g, area, y, a, b, kind in jobs:
        fetch(client, out, kind, g, area, y, a, b, failed)

    print("Finished. Failed requests: %d" % len(failed))
    for f in failed:
        print("   ", f)
    if failed:
        print("Run the same command again to retry them. If every request fails, check ~/.cdsapirc "
              "and that you accepted the dataset terms on the CDS website.")
        sys.exit(2)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopped by you. Run the same command again to continue.")
        sys.exit(1)
