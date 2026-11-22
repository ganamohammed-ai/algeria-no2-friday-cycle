#!/usr/bin/env python3
"""Download TROPOMI (Sentinel-5P) Level-2 NO2 orbit files from the Copernicus Data
Space Ecosystem (CDSE), cut a small box around each study site, keep the cut
pieces, and delete the big original file.

Typical use (see README_AR.txt for the Arabic step-by-step guide):

    python tropomi_download.py --count-only                 # how many files / GB?
    python tropomi_download.py --inspect                    # download ONE file, show its variables
    python tropomi_download.py --start 2024-01-01 --end 2024-01-31   # one-month test
    python tropomi_download.py                              # everything 2019-08-07 .. 2025-12-31

The script is resumable: run it again and it skips orbits already done.
Your CDSE password is asked at run time (or taken from the environment
variables CDSE_USER / CDSE_PASS) and is never written to disk.
"""
import argparse
import csv
import datetime as dt
import getpass
import json
import os
import queue
import re
import sys
import threading
import time
from pathlib import Path

import numpy as np
import requests

from sites import SITES, parse_site_list, site_box, HALF_BOX_DEG

CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URLS = [
    "https://download.dataspace.copernicus.eu/odata/v1/Products({id})/$value",
    "https://zipper.dataspace.copernicus.eu/odata/v1/Products({id})/$value",
]
TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"

# S5P_OFFL_L2__NO2____20240101T115000_20240101T133130_32000_03_020600_20240103T...nc
NAME_RE = re.compile(
    r"S5P_(OFFL|RPRO)_L2__NO2____(\d{8}T\d{6})_(\d{8}T\d{6})_(\d{5})_(\d{2})_(\d{6})_(\d{8}T\d{6})(?:\.nc)?$"
)

# Variables to keep. Required ones must exist; the others are saved when found.
REQUIRED = ["latitude", "longitude", "nitrogendioxide_tropospheric_column", "qa_value"]
OPTIONAL_PIXEL = [
    "nitrogendioxide_tropospheric_column_precision",
    "air_mass_factor_troposphere",
    "air_mass_factor_total",
    "nitrogendioxide_stratospheric_column",
    "nitrogendioxide_total_column",
    "cloud_fraction_crb",
    "cloud_radiance_fraction_nitrogendioxide_window",
    "surface_albedo_nitrogendioxide_window",
    "surface_pressure",
    "solar_zenith_angle",
    "viewing_zenith_angle",
    "solar_azimuth_angle",
    "viewing_azimuth_angle",
    "tm5_tropopause_layer_index",
]
CORNER_VARS = ["latitude_bounds", "longitude_bounds"]
AK_VAR = "averaging_kernel"
LAYER_VARS = ["tm5_constant_a", "tm5_constant_b"]  # (layer, 2), not per pixel


# ----------------------------------------------------------------------------
# Catalogue search
# ----------------------------------------------------------------------------
class Auth:
    """Keeps a valid CDSE access token (tokens last about 10 minutes)."""

    def __init__(self, user, password):
        self.user, self.password = user, password
        self.token, self.t0 = None, 0.0
        self.lock = threading.Lock()

    def get(self):
        with self.lock:
            return self._get()

    def _get(self):
        if self.token is None or time.time() - self.t0 > 420:
            r = requests.post(
                TOKEN_URL,
                data={
                    "client_id": "cdse-public",
                    "grant_type": "password",
                    "username": self.user,
                    "password": self.password,
                },
                timeout=60,
            )
            if r.status_code != 200:
                raise SystemExit(
                    "Login to CDSE failed (HTTP %s). Check user name / password, "
                    "and that you can log in at https://dataspace.copernicus.eu" % r.status_code
                )
            self.token, self.t0 = r.json()["access_token"], time.time()
        return self.token


def polygon_wkt(box):
    x0, y0, x1, y1 = box
    return "POLYGON((%.4f %.4f,%.4f %.4f,%.4f %.4f,%.4f %.4f,%.4f %.4f))" % (
        x0, y0, x1, y0, x1, y1, x0, y1, x0, y0)


def odata_filter(box, start, end):
    return (
        "Collection/Name eq 'SENTINEL-5P' and "
        "Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and "
        "att/OData.CSC.StringAttribute/Value eq 'L2__NO2___') and "
        "OData.CSC.Intersects(area=geography'SRID=4326;%s') and "
        "ContentDate/Start gt %sT00:00:00.000Z and ContentDate/Start lt %sT00:00:00.000Z"
        % (polygon_wkt(box), start.isoformat(), end.isoformat())
    )


def search(box, start, end):
    """All L2 NO2 products intersecting `box` with start <= date < end."""
    out, skip, top = [], 0, 200
    while True:
        params = {"$filter": odata_filter(box, start, end), "$orderby": "ContentDate/Start asc",
                  "$top": str(top), "$skip": str(skip)}
        for attempt in range(4):
            try:
                r = requests.get(CATALOGUE, params=params, timeout=120)
                r.raise_for_status()
                break
            except requests.RequestException as e:
                if attempt == 3:
                    raise
                time.sleep(5 * (attempt + 1))
        batch = r.json().get("value", [])
        out.extend(batch)
        if len(batch) < top:
            return out
        skip += top


def month_ranges(start, end):
    """Yield (first_day, next_first_day) chunks of one month, clipped to [start, end+1)."""
    cur = start
    stop = end + dt.timedelta(days=1)
    while cur < stop:
        nxt = (cur.replace(day=1) + dt.timedelta(days=32)).replace(day=1)
        yield cur, min(nxt, stop)
        cur = nxt


def _rings(geo):
    """Exterior rings [(x, y), ...] from a GeoJSON Polygon / MultiPolygon."""
    if isinstance(geo, str):
        try:
            geo = json.loads(geo)
        except ValueError:
            return []
    if not isinstance(geo, dict):
        return []
    t, c = geo.get("type"), geo.get("coordinates")
    if not c:
        return []
    if t == "Polygon":
        return [[(p[0], p[1]) for p in c[0]]]
    if t == "MultiPolygon":
        return [[(p[0], p[1]) for p in poly[0]] for poly in c]
    return []


def _clip_to_box(poly, x0, y0, x1, y1):
    """Sutherland-Hodgman clipping of a polygon to a rectangle."""
    def clip_edge(pts, inside, inter):
        out = []
        if not pts:
            return out
        prev = pts[-1]
        for cur in pts:
            if inside(cur):
                if not inside(prev):
                    out.append(inter(prev, cur))
                out.append(cur)
            elif inside(prev):
                out.append(inter(prev, cur))
            prev = cur
        return out

    def at_x(xc):
        return lambda a, b: (xc, a[1] + (b[1] - a[1]) * (xc - a[0]) / (b[0] - a[0]))

    def at_y(yc):
        return lambda a, b: (a[0] + (b[0] - a[0]) * (yc - a[1]) / (b[1] - a[1]), yc)

    pts = clip_edge(poly, lambda p: p[0] >= x0, at_x(x0))
    pts = clip_edge(pts, lambda p: p[0] <= x1, at_x(x1))
    pts = clip_edge(pts, lambda p: p[1] >= y0, at_y(y0))
    pts = clip_edge(pts, lambda p: p[1] <= y1, at_y(y1))
    return pts


def _area(pts):
    a = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def _point_in_ring(pt, ring):
    x, y = pt
    inside = False
    n = len(ring)
    for i in range(n):
        xa, ya = ring[i]
        xb, yb = ring[(i + 1) % n]
        if (ya > y) != (yb > y) and x < (xb - xa) * (y - ya) / (yb - ya) + xa:
            inside = not inside
    return inside


def footprint_overlap(geo, box, point):
    """(fraction of the box covered by the orbit footprint, site point inside it)
    or None when the catalogue gave no footprint."""
    rings = _rings(geo)
    if not rings:
        return None
    x0, y0, x1, y1 = box
    area = sum(_area(_clip_to_box(r, x0, y0, x1, y1)) for r in rings)
    frac = min(area / ((x1 - x0) * (y1 - y0)), 1.0)
    return frac, any(_point_in_ring(point, r) for r in rings)


def build_catalogue(site_names, start, end, half, min_overlap=0.0):
    """(orbit, date) -> info dict with the list of sites it usefully covers.
    A product is used for a site only if its footprint covers at least `min_overlap`
    of the site box AND contains the site centre (skipped when the footprint is unknown)."""
    cat, skipped = {}, 0
    for s in site_names:
        box = site_box(s, half)
        lat, lon = SITES[s]
        n_site = n_kept = 0
        for a, b in month_ranges(start, end):
            for p in search(box, a, b):
                m = NAME_RE.match(p["Name"][:-3] if p["Name"].endswith(".nc") else p["Name"])
                if not m:
                    continue
                n_site += 1
                ov = footprint_overlap(p.get("GeoFootprint"), box, (lon, lat))
                if min_overlap > 0 and ov is not None and (not ov[1] or ov[0] < min_overlap):
                    skipped += 1
                    continue
                n_kept += 1
                mode, t0, t1, orbit, coll, proc, prod = m.groups()
                key = (orbit, t0[:8])
                info = cat.get(key)
                if info is None or (mode == "RPRO" and info["mode"] != "RPRO"):
                    sites_prev = info["sites"] if info else set()
                    info = {"id": p["Id"], "name": p["Name"], "mode": mode, "start": t0, "end": t1,
                            "orbit": orbit, "processor": proc, "size": p.get("ContentLength", 0) or 0,
                            "sites": sites_prev}
                    cat[key] = info
                info["sites"].add(s)
        print("  searched %-15s %d products, %d useful" % (s, n_site, n_kept), flush=True)
    return cat, skipped


# ----------------------------------------------------------------------------
# Download
# ----------------------------------------------------------------------------
def download(auth, info, dest):
    part = Path(str(dest) + ".part")
    last_err = None
    for url_tpl in DOWNLOAD_URLS:
        url = url_tpl.format(id=info["id"])
        for attempt in range(3):
            try:
                headers = {"Authorization": "Bearer " + auth.get()}
                with requests.get(url, headers=headers, stream=True, timeout=(30, 300)) as r:
                    if r.status_code in (401, 403):
                        auth.token = None
                        raise requests.RequestException("HTTP %s" % r.status_code)
                    r.raise_for_status()
                    with open(part, "wb") as f:
                        for chunk in r.iter_content(chunk_size=1 << 20):
                            f.write(chunk)
                if info["size"] and part.stat().st_size != int(info["size"]):
                    raise requests.RequestException("size mismatch")
                part.replace(dest)
                return True
            except requests.RequestException as e:
                last_err = e
                time.sleep(10 * (attempt + 1))
    print("    download failed: %s" % last_err, flush=True)
    if part.exists():
        part.unlink()
    return False


# ----------------------------------------------------------------------------
# Reading and cutting a TROPOMI file
# ----------------------------------------------------------------------------
def walk_vars(group, prefix=""):
    for name, var in group.variables.items():
        yield prefix + name, name, var
    for gname, g in group.groups.items():
        yield from walk_vars(g, prefix + gname + "/")


def find_var(ds, name):
    for _full, short, var in walk_vars(ds):
        if short == name:
            return var
    return None


def to_array(var):
    """Variable -> numpy array without the leading size-1 'time' axis.
    Floats: masked -> NaN. Integers: masked -> -1."""
    a = var[:]
    if np.ma.isMaskedArray(a):
        if a.dtype.kind == "f":
            a = a.astype("float32").filled(np.nan)
        elif a.dtype.kind in "iu":
            a = a.filled(-1)
        else:
            a = a.filled("")
    if var.dimensions and var.dimensions[0] == "time" and a.shape[0] == 1:
        a = a[0]
    return a


def print_tree(ds):
    for full, _short, var in walk_vars(ds):
        print("  %-90s %s %s" % (full, var.dtype, var.shape))


def cut_file(path, site_names, outdir, half, keep_ak, info):
    """Cut a box around each site and save .npz files. Returns {site: n_valid_pixels}."""
    from netCDF4 import Dataset  # imported here so --count-only works without it

    saved = {}
    with Dataset(path) as ds:
        for r in REQUIRED:
            if find_var(ds, r) is None:
                raise RuntimeError("required variable not found: " + r)
        lat = to_array(find_var(ds, "latitude"))
        lon = to_array(find_var(ds, "longitude"))
        no2 = to_array(find_var(ds, "nitrogendioxide_tropospheric_column"))
        qa = to_array(find_var(ds, "qa_value"))
        tvar = find_var(ds, "time_utc")
        times = np.asarray(to_array(tvar)).astype(str) if tvar is not None else None

        cache = {}  # variables read lazily, once

        def pix(name):
            if name not in cache:
                v = find_var(ds, name)
                cache[name] = None if v is None else to_array(v)
            return cache[name]

        date = info["start"][:8]
        for s in site_names:
            x0, y0, x1, y1 = site_box(s, half)
            m = (lat >= y0) & (lat <= y1) & (lon >= x0) & (lon <= x1) & np.isfinite(no2)
            if not m.any():
                continue
            sc, gp = np.nonzero(m)
            d = {
                "scanline": sc.astype("int32"),
                "ground_pixel": gp.astype("int16"),
                "latitude": lat[sc, gp].astype("float32"),
                "longitude": lon[sc, gp].astype("float32"),
                "no2_trop_mol_m2": no2[sc, gp].astype("float32"),
                "qa_value": qa[sc, gp].astype("float32"),
            }
            if times is not None:
                d["time_utc"] = times[sc]
            for name in OPTIONAL_PIXEL:
                a = pix(name)
                if a is not None and a.ndim == 2:
                    d[name] = a[sc, gp]
            for name in CORNER_VARS:
                a = pix(name)
                if a is not None and a.ndim == 3:
                    d[name] = a[sc, gp, :].astype("float32")
            if keep_ak:
                a = pix(AK_VAR)
                if a is not None and a.ndim == 3:
                    d[AK_VAR] = a[sc, gp, :].astype("float32")
                for name in LAYER_VARS:
                    a = pix(name)
                    if a is not None:
                        d[name] = np.asarray(a)
            d["product_name"] = np.array(info["name"])
            d["processing_mode"] = np.array(info["mode"])
            d["processor_version"] = np.array(info["processor"])
            d["orbit"] = np.array(info["orbit"])
            sdir = Path(outdir) / s
            sdir.mkdir(parents=True, exist_ok=True)
            target = sdir / ("%s_%s.npz" % (date, info["orbit"]))
            tmp = sdir / ("%s_%s.tmp.npz" % (date, info["orbit"]))
            np.savez_compressed(tmp, **d)
            os.replace(tmp, target)
            saved[s] = int(((qa[sc, gp] > 0.75) & np.isfinite(no2[sc, gp])).sum())
    return saved


# ----------------------------------------------------------------------------
# Progress file (so a re-run skips finished orbits)
# ----------------------------------------------------------------------------
def load_done(path):
    done = set()
    if path.exists():
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row["status"] in ("done", "nodata"):
                    done.add(row["name"])
    return done


def log_progress(path, info, status, saved):
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["name", "orbit", "date", "mode", "processor", "status", "sites_valid_pixels"])
        w.writerow([info["name"], info["orbit"], info["start"][:8], info["mode"], info["processor"],
                    status, ";".join("%s=%d" % kv for kv in sorted(saved.items()))])


# ----------------------------------------------------------------------------
def get_auth():
    user = os.environ.get("CDSE_USER") or input("CDSE e-mail/user name: ").strip()
    pw = os.environ.get("CDSE_PASS") or getpass.getpass("CDSE password (hidden, not saved): ")
    return Auth(user, pw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="2019-08-07", help="first day, YYYY-MM-DD (default 2019-08-07)")
    ap.add_argument("--end", default="2025-12-31", help="last day, YYYY-MM-DD (default 2025-12-31)")
    ap.add_argument("--sites", default="default", help="comma list, e.g. algiers,constantine (default: all)")
    ap.add_argument("--half-box", type=float, default=HALF_BOX_DEG, help="box radius in degrees of latitude")
    ap.add_argument("--min-overlap", type=float, default=0.4,
                    help="skip orbit files whose footprint covers less than this fraction of a site box "
                         "or misses the site centre (default 0.4; 0 = keep everything)")
    ap.add_argument("--workers", type=int, default=3,
                    help="parallel downloads (default 3; the server allows at most 4)")
    ap.add_argument("--outdir", default="tropomi_clips", help="folder for the cut pieces")
    ap.add_argument("--rawdir", default="tropomi_raw_tmp", help="temporary folder for the big files")
    ap.add_argument("--no-ak", action="store_true", help="do not keep averaging kernels (smaller output)")
    ap.add_argument("--keep-raw", action="store_true", help="do not delete the big original files")
    ap.add_argument("--count-only", action="store_true", help="only count files and size, download nothing")
    ap.add_argument("--inspect", action="store_true", help="download ONE file and print its variables")
    args = ap.parse_args()

    start = dt.date.fromisoformat(args.start)
    end = dt.date.fromisoformat(args.end)
    sites = parse_site_list(args.sites)
    workers = max(1, min(4, args.workers))
    print("Sites: %s | %s .. %s | box radius %.2f deg | min overlap %.2f | workers %d" % (
        ", ".join(sites), start, end, args.half_box, args.min_overlap, workers))

    print("Searching the catalogue ...", flush=True)
    cat, skipped = build_catalogue(sites, start, end, args.half_box, args.min_overlap)
    items = sorted(cat.values(), key=lambda i: (i["start"], i["orbit"]))
    total_gb = sum(int(i["size"]) for i in items) / 1e9
    print("Skipped by the footprint filter (edge-only or missing the site): %d site-products" % skipped)
    print("Orbit files to use: %d   total size to transfer: %.0f GB" % (len(items), total_gb))
    by_mode = {}
    for i in items:
        by_mode[(i["mode"], i["processor"])] = by_mode.get((i["mode"], i["processor"]), 0) + 1
    for (mode, proc), n in sorted(by_mode.items()):
        print("   %s processor %s : %d files" % (mode, proc, n))
    if args.count_only:
        return

    outdir, rawdir = Path(args.outdir), Path(args.rawdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rawdir.mkdir(parents=True, exist_ok=True)
    progress = outdir / "progress.csv"
    done = load_done(progress)
    auth = get_auth()
    auth.get()  # fail early if the login is wrong

    todo = [i for i in items if i["name"] not in done]
    print("To do now: %d (already done: %d)" % (len(todo), len(items) - len(todo)), flush=True)

    def raw_path(info):
        return rawdir / (info["name"] if info["name"].endswith(".nc") else info["name"] + ".nc")

    if args.inspect:
        info = todo[0]
        raw = raw_path(info)
        if not raw.exists() and not download(auth, info, raw):
            raise SystemExit("download failed")
        from netCDF4 import Dataset
        with Dataset(raw) as ds:
            print_tree(ds)
        print("\nRaw file kept at: %s" % raw)
        return

    q = queue.Queue()
    for info in todo:
        q.put(info)
    st = {"done": 0, "fail": 0, "consec": 0, "stop": False, "abort": False}
    lock, cut_lock = threading.Lock(), threading.Lock()
    t_begin, total = time.time(), len(todo)

    def eta_text():
        k = st["done"]
        if k < 3:
            return ""
        return "  ETA %.1f h" % ((time.time() - t_begin) / k * (total - k) / 3600)

    def work():
        while not st["stop"]:
            try:
                info = q.get_nowait()
            except queue.Empty:
                return
            raw = raw_path(info)
            if not raw.exists() and not download(auth, info, raw):
                with lock:
                    st["fail"] += 1
                    st["done"] += 1
                    log_progress(progress, info, "download_failed", {})
                    print("[%d/%d] %s orbit %s  DOWNLOAD FAILED (will retry on re-run)" % (
                        st["done"], total, info["start"][:8], info["orbit"]), flush=True)
                continue
            try:
                with cut_lock:  # netCDF4 is not thread-safe: cut one file at a time
                    saved = cut_file(raw, sorted(info["sites"]), outdir, args.half_box, not args.no_ak, info)
                with lock:
                    st["consec"] = 0
                    st["done"] += 1
                    log_progress(progress, info, "done" if saved else "nodata", saved)
                    print("[%d/%d] %s %s orbit %s saved: %s%s" % (
                        st["done"], total, info["start"][:8], info["mode"], info["orbit"],
                        ", ".join("%s(%d valid)" % kv for kv in saved.items()) or "nothing in boxes",
                        eta_text()), flush=True)
                if not args.keep_raw and raw.exists():
                    raw.unlink()
            except Exception as e:  # the orbit is retried on the next run
                with lock:
                    st["fail"] += 1
                    st["done"] += 1
                    st["consec"] += 1
                    log_progress(progress, info, "cut_failed: %s" % e, {})
                    print("[%d/%d] %s orbit %s  CUTTING FAILED: %s" % (
                        st["done"], total, info["start"][:8], info["orbit"], e), flush=True)
                    if st["consec"] >= 3 and not st["abort"]:
                        st["abort"] = st["stop"] = True
                        print("\nSTOPPING: 3 files in a row could not be read/cut. The file layout is probably "
                              "not what this script expects. The failing files are kept in %s . "
                              "Send the error text above to Claude." % rawdir, flush=True)

    threads = [threading.Thread(target=work, daemon=True) for _ in range(workers)]
    for t in threads:
        t.start()
        time.sleep(2)
    while any(t.is_alive() for t in threads):
        time.sleep(1)

    if st["abort"]:
        sys.exit(3)
    print("Finished. Failed orbits: %d (run the script again to retry them)." % st["fail"])
    sys.exit(2 if st["fail"] else 0)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopped by you. Run the same command again to continue where it stopped.")
        sys.exit(1)
