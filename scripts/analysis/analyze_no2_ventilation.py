#!/usr/bin/env python3
"""NO2 (TROPOMI, from Google Earth Engine) versus atmospheric ventilation (ERA5), six Algerian sites.

Three stages, each can be run on its own:

  1) tables : read the GEE GeoTIFFs (one band per day, NO2 in umol/m2, -9999 = no data) and build a
              daily table per site: mean NO2 within 10 km of the site ("inner"), mean in the 50-100 km
              ring ("ring"), box mean and the number of valid pixels.
  2) merge  : read ERA5 (the CSV files of ERA5_GEE_export.js, or the NetCDF files of era5_download.py;
              12:00 and 13:00 UTC) and add
              wind, boundary-layer height (BLH), ventilation coefficient and temperature to each day.
  3) model  : per-site regression of ln(NO2) on ln(ventilation) with controls, plus a binned plot.

Examples (run from the folder that holds the GeoTIFFs):
    python analyze_no2_ventilation.py --stage tables --tif-dir . --out out
    python analyze_no2_ventilation.py --stage all    --tif-dir . --era5-dir era5 --out out
    python analyze_no2_ventilation.py --stage tables --tif-dir . --out out --sites algiers --years 2022

Packages: numpy, tifffile and imagecodecs (the GeoTIFFs are LZW-compressed), netCDF4 (for ERA5), matplotlib (plots).
    pip install numpy tifffile imagecodecs netCDF4 matplotlib

Definitions fixed here (change them at the top if the analysis plan says otherwise):
  * date = UTC date of the GEE band; ERA5 = mean of 12:00 and 13:00 UTC (the satellite passes at ~12:30 UTC).
  * ventilation coefficient VC = BLH x wind speed. Two versions:
      VC10    = BLH x |10 m wind|
      VCLAYER = BLH x |mean wind vector inside the boundary layer| (pressure levels between the surface and
                the BLH top; the 10 m wind counts as the lowest level). BLH top pressure =
                surface pressure x exp(-g x BLH / (R x T2m)).
  * a day is used only if at least MIN_FRAC of the pixels of the inner disc AND of the ring are valid.
"""
import argparse
import csv
import datetime as dt
import math
import os
import re
import sys

import numpy as np

# ------------------------------- SETTINGS -----------------------------------
SITES = {  # (lat, lon): same as sites.py
    "algiers": (36.75, 3.06),
    "oran": (35.70, -0.63),
    "arzew": (35.85, -0.31),
    "constantine": (36.37, 6.61),
    "setif": (36.19, 5.41),
    "hassi_messaoud": (31.68, 6.07),
}
# Control (background) sites, chosen from static layers only (select_control.py --amended). Not part of the default site list.
CONTROL_SITES = {"control_tell": (35.0, 5.25), "control_highland": (34.5, -1.0)}
DEFAULT_SITES = list(SITES)
SITES.update(CONTROL_SITES)
HALF_DEG = 1.5            # box radius used for the GEE export (degrees of latitude)
INNER_KM = 10.0
RING_KM = (50.0, 100.0)
MIN_FRAC = 0.30           # minimum share of valid pixels in the inner disc and in the ring
NODATA = -9999.0
FIRST_DAY = dt.date(2019, 8, 7)    # same as the GEE script
LAST_DAY = dt.date(2025, 12, 31)
YEARS = list(range(2019, 2026))
G = 9.80665
R_AIR = 287.05
NW_LAG = 7                # Newey-West lag (days)
# -----------------------------------------------------------------------------


def year_days(y):
    a = max(dt.date(y, 1, 1), FIRST_DAY)
    b = min(dt.date(y, 12, 31), LAST_DAY)
    n = (b - a).days + 1
    return [a + dt.timedelta(days=i) for i in range(n)]


# ================================ stage 1 ====================================
def read_geotiff(path):
    """Return (array H x W x B float32 with NaN for no data, lon0, lat0, dx, dy) with lon0/lat0 = top-left corner."""
    try:
        import tifffile
    except ImportError:
        raise SystemExit("Missing package. Run:  pip install tifffile imagecodecs")
    with tifffile.TiffFile(path) as f:
        page = f.pages[0]
        planar_separate = (str(page.planarconfig).endswith("SEPARATE") or page.planarconfig == 2)
        scale = page.tags["ModelPixelScaleTag"].value
        tie = page.tags["ModelTiepointTag"].value
        try:
            arr = f.asarray()
        except ValueError as e:
            raise SystemExit("Cannot decode %s (%s). Run:  pip install imagecodecs" % (path, e))
    dx, dy = float(scale[0]), float(scale[1])
    lon0 = float(tie[3]) - float(tie[0]) * dx
    lat0 = float(tie[4]) + float(tie[1]) * dy
    if arr.ndim == 2:
        arr = arr[:, :, None]
    elif planar_separate:
        arr = np.moveaxis(arr, 0, -1)       # (bands, H, W) -> (H, W, bands)
    arr = arr.astype("float32", copy=False)
    arr = np.where(arr <= -9000, np.nan, arr)
    return arr, lon0, lat0, dx, dy


def disc_masks(shape_hw, lon0, lat0, dx, dy, site):
    lat_s, lon_s = SITES[site]
    H, W = shape_hw
    lon = lon0 + (np.arange(W) + 0.5) * dx
    lat = lat0 - (np.arange(H) + 0.5) * dy
    LON, LAT = np.meshgrid(lon, lat)
    ykm = (LAT - lat_s) * 110.574
    xkm = (LON - lon_s) * 111.320 * math.cos(math.radians(lat_s))
    r = np.hypot(xkm, ykm)
    return r <= INNER_KM, (r >= RING_KM[0]) & (r <= RING_KM[1]), np.ones_like(r, bool)


def masked_mean(arr, mask):
    sub = arr[mask]                       # (npix, B)
    n = np.sum(~np.isnan(sub), axis=0)
    with np.errstate(invalid="ignore"):
        s = np.nansum(sub, axis=0)
        mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return mean, n, mask.sum()


def build_no2_rows(tif_path, site, year):
    arr, lon0, lat0, dx, dy = read_geotiff(tif_path)
    days = year_days(year)
    if arr.shape[2] != len(days):
        raise SystemExit("%s has %d bands but %d days were expected for %d" % (tif_path, arr.shape[2], len(days), year))
    inner, ring, box = disc_masks(arr.shape[:2], lon0, lat0, dx, dy, site)
    mi, ni, ti = masked_mean(arr, inner)
    mr, nr, tr = masked_mean(arr, ring)
    mb, nb, tb = masked_mean(arr, box)
    rows = []
    for k, d in enumerate(days):
        ok = (ni[k] >= MIN_FRAC * ti) and (nr[k] >= MIN_FRAC * tr)
        rows.append({
            "site": site, "date": d.isoformat(),
            "no2_inner": mi[k], "inner_n": int(ni[k]), "inner_total": int(ti),
            "no2_ring": mr[k], "ring_n": int(nr[k]), "ring_total": int(tr),
            "no2_box": mb[k], "box_n": int(nb[k]),
            "no2_delta": (mi[k] - mr[k]) if ok else np.nan,
            "usable": int(ok),
        })
    return rows


def stage_tables(tif_dir, sites, years, out):
    os.makedirs(out, exist_ok=True)
    allrows = []
    for site in sites:
        for y in years:
            p = os.path.join(tif_dir, "NO2_daily_%s_%d.tif" % (site, y))
            if not os.path.exists(p):
                print("  missing file:", p)
                continue
            rows = build_no2_rows(p, site, y)
            u = sum(r["usable"] for r in rows)
            print("  %-15s %d: %3d days, %3d usable (%.0f%%)" % (site, y, len(rows), u, 100.0 * u / len(rows)), flush=True)
            allrows += rows
    path = os.path.join(out, "no2_daily_table.csv")
    write_csv(path, allrows)
    print("wrote", path, "(%d rows)" % len(allrows))
    return allrows


def write_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and math.isnan(v)) else v) for k, v in r.items()})


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            if k in ("site", "date"):
                continue
            r[k] = float(v) if v != "" else float("nan")
    return rows


# ================================ stage 2 ====================================
class NC:
    """Minimal NetCDF reader: netCDF4 if present, otherwise scipy (classic files only; used in tests)."""

    def __init__(self, path):
        self.path = path
        try:
            import netCDF4
            self.ds = netCDF4.Dataset(path)
            self.ds.set_auto_maskandscale(False)
            self.kind = "nc4"
        except ImportError:
            from scipy.io import netcdf_file
            self.ds = netcdf_file(path, mmap=False)
            self.kind = "scipy"

    def names(self):
        return list(self.ds.variables.keys())

    def get(self, name):
        v = self.ds.variables[name]
        a = np.array(v[:])
        attrs = {k: getattr(v, k) for k in ("scale_factor", "add_offset", "_FillValue", "missing_value", "units") if hasattr(v, k)}
        if "_FillValue" in attrs and a.dtype.kind in "iuf":
            a = np.where(a == attrs["_FillValue"], np.nan, a.astype("float64"))
        if "scale_factor" in attrs:
            a = a.astype("float64") * float(np.squeeze(attrs["scale_factor"])) + float(np.squeeze(attrs.get("add_offset", 0.0)))
        return a, attrs

    def first(self, *cands):
        for c in cands:
            if c in self.ds.variables:
                return c
        return None

    def close(self):
        try:
            self.ds.close()
        except Exception:
            pass


def nc_times(nc):
    tname = nc.first("valid_time", "time")
    a, at = nc.get(tname)
    units = at.get("units", b"")
    units = units.decode() if isinstance(units, bytes) else str(units)
    m = re.match(r"(\w+) since (\d{4}-\d{2}-\d{2})(?:[ T](\d{2}):(\d{2}):(\d{2}))?", units)
    if not m:
        raise SystemExit("Cannot parse time units %r in %s" % (units, nc.path))
    unit, d0 = m.group(1), dt.datetime.strptime(m.group(2), "%Y-%m-%d")
    mult = {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}[unit]
    return [d0 + dt.timedelta(seconds=float(x) * mult) for x in np.ravel(a)]


def bilinear(grid_lat, grid_lon, field, lat, lon):
    """field[..., nlat, nlon]; lat/lon vectors may be in any order. Returns field interpolated at the point."""
    la = np.asarray(grid_lat, float)
    lo = np.asarray(grid_lon, float)
    if la[0] > la[-1]:
        la = la[::-1]
        field = field[..., ::-1, :]
    i = int(np.clip(np.searchsorted(la, lat) - 1, 0, len(la) - 2))
    j = int(np.clip(np.searchsorted(lo, lon) - 1, 0, len(lo) - 2))
    ty = (lat - la[i]) / (la[i + 1] - la[i])
    tx = (lon - lo[j]) / (lo[j + 1] - lo[j])
    f00, f01 = field[..., i, j], field[..., i, j + 1]
    f10, f11 = field[..., i + 1, j], field[..., i + 1, j + 1]
    return (1 - ty) * ((1 - tx) * f00 + tx * f01) + ty * ((1 - tx) * f10 + tx * f11)


def squeeze_extra(a, expect_ndim):
    """CDS files for recent months can carry an 'expver' axis; collapse it."""
    while a.ndim > expect_ndim:
        a = np.nanmax(a, axis=1) if a.shape[1] <= 5 else a[:, 0]
    return a


def nc_members(path, tmp):
    """CDS sometimes delivers a ZIP holding several NetCDF files (e.g. instant and accumulated variables).
    Return the list of NetCDF paths to read: the file itself, or the members extracted into tmp."""
    import tempfile
    import zipfile
    if not zipfile.is_zipfile(path):
        return [path]
    d = tempfile.mkdtemp(dir=tmp)
    with zipfile.ZipFile(path) as z:
        z.extractall(d)
        return [os.path.join(d, n) for n in sorted(z.namelist()) if n.lower().endswith(".nc")]


def read_era5_set(sl_paths, pl_paths, sites):
    """Return {(site, date): dict(variables)} from any set of ERA5 files (yearly, quarterly, monthly).
    Single-level and pressure-level files are matched by time stamp, so they need not be cut the same way.
    The 12 and 13 UTC fields are averaged; the accumulated fields (ssrd, tp) are the 12->13 UTC increment."""
    import shutil
    import tempfile
    out = {}
    tmp = tempfile.mkdtemp(prefix="era5_unzip_")
    try:
        return _read_era5_set(sl_paths, pl_paths, sites, tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _read_era5_set(sl_paths, pl_paths, sites, tmp):
    out = {}
    keys = (("u10", ("u10",)), ("v10", ("v10",)), ("blh", ("blh",)), ("t2m", ("t2m",)),
            ("sp", ("sp",)), ("ssrd", ("ssrd",)), ("tp", ("tp",)))
    t, lat, lon, parts = [], None, None, {}
    for p in sl_paths:
        got = {}                                   # key -> (times, array) for this file
        for mp in nc_members(p, tmp):
            f = NC(mp)
            tm = nc_times(f)
            if lat is None and f.first("latitude", "lat"):
                lat = f.get(f.first("latitude", "lat"))[0]
                lon = f.get(f.first("longitude", "lon"))[0]
            for key, cands in keys:
                nm = f.first(*cands)
                if nm is not None:
                    got[key] = (tm, squeeze_extra(f.get(nm)[0], 3))   # (time, lat, lon)
            f.close()
        if "u10" not in got:
            print("  skipped (no u10 found):", p)
            continue
        tfile = got["u10"][0]
        for key, (tm, arr) in got.items():
            if len(tm) != len(tfile):
                print("  dropped %s from %s (different time axis)" % (key, os.path.basename(p)))
                continue
            parts.setdefault(key, []).append(arr)
        t += tfile
    if not t:
        return out
    var = {k: np.concatenate(v, axis=0) for k, v in parts.items()}
    pl = None
    tpl, Ul, Vl = [], [], []
    plev = plat = plon = None
    for p in pl_paths:
        if not os.path.exists(p):
            continue
        for mp in nc_members(p, tmp):
            f = NC(mp)
            if f.first("u") is None:
                f.close()
                continue
            tpl += nc_times(f)
            if plev is None:
                plev = np.ravel(f.get(f.first("pressure_level", "level"))[0]).astype(float)   # hPa
                plat = f.get(f.first("latitude", "lat"))[0]
                plon = f.get(f.first("longitude", "lon"))[0]
            Ul.append(squeeze_extra(f.get(f.first("u"))[0], 4))        # (time, level, lat, lon)
            Vl.append(squeeze_extra(f.get(f.first("v"))[0], 4))
            f.close()
    if tpl:
        pl = True
        U = np.concatenate(Ul, axis=0)
        V = np.concatenate(Vl, axis=0)
        tmap = {x: i for i, x in enumerate(tpl)}
    days = sorted({x.date() for x in t})
    nan = float("nan")
    for site in sites:
        la, lo = SITES[site]
        s = {k: np.array([bilinear(lat, lon, v[i], la, lo) for i in range(v.shape[0])]) for k, v in var.items()}
        if pl:
            su = np.array([bilinear(plat, plon, U[i], la, lo) for i in range(U.shape[0])])   # (time, level)
            sv = np.array([bilinear(plat, plon, V[i], la, lo) for i in range(V.shape[0])])
        for d in days:
            idx = [i for i, x in enumerate(t) if x.date() == d and x.hour in (12, 13)]
            if not idx:
                continue
            i12 = [i for i in idx if t[i].hour == 12]
            i13 = [i for i in idx if t[i].hour == 13]
            m = {k: float(np.nanmean(v[idx])) for k, v in s.items() if k not in ("ssrd", "tp")}
            rec = {}
            rec["u10"], rec["v10"] = m["u10"], m["v10"]
            rec["wspd10"] = math.hypot(m["u10"], m["v10"])
            rec["wdir10"] = (math.degrees(math.atan2(-m["u10"], -m["v10"])) + 360.0) % 360.0   # direction wind comes FROM
            rec["blh"] = m["blh"]
            rec["t2m"] = m["t2m"] - 273.15
            rec["sp_hpa"] = m["sp"] / 100.0
            # accumulated since the forecast start: the hourly amount is the 13 UTC minus the 12 UTC value
            rec["ssrd_jm2"] = float(s["ssrd"][i13[0]] - s["ssrd"][i12[0]]) if ("ssrd" in s and i12 and i13) else nan
            rec["tp_mm"] = float(s["tp"][i13[0]] - s["tp"][i12[0]]) * 1000.0 if ("tp" in s and i12 and i13) else nan
            rec["vc10"] = rec["blh"] * rec["wspd10"]
            for k in ("wspd_layer", "wdir_layer", "n_levels_in_blh", "vclayer"):
                rec[k] = nan
            pidx = [tmap[t[i]] for i in idx if t[i] in tmap] if pl else []
            if pidx:
                ptop = m["sp"] / 100.0 * math.exp(-G * m["blh"] / (R_AIR * m["t2m"]))
                uu = [m["u10"]]
                vv = [m["v10"]]
                for li, p in enumerate(plev):
                    if ptop <= p <= m["sp"] / 100.0:
                        uu.append(float(np.nanmean(su[pidx, li])))
                        vv.append(float(np.nanmean(sv[pidx, li])))
                ml, mv = float(np.mean(uu)), float(np.mean(vv))
                rec["wspd_layer"] = math.hypot(ml, mv)
                rec["wdir_layer"] = (math.degrees(math.atan2(-ml, -mv)) + 360.0) % 360.0
                rec["n_levels_in_blh"] = len(uu) - 1
                rec["vclayer"] = rec["blh"] * rec["wspd_layer"]
            out[(site, d.isoformat())] = rec
    return out


def profile_layer_wind(u10, v10, u100, v100, u850, v850, u500, v500, sp_pa, t2m_k, blh):
    """Mean wind vector from the ground up to the boundary-layer top, from a piecewise-linear profile through
    the 10 m, 100 m, 850 hPa and 500 hPa winds (heights above ground from the hypsometric equation).
    Approximation used when the full set of pressure levels is not available."""
    nodes = [(10.0, u10, v10), (100.0, u100, v100)]
    for p, uu, vv in ((850.0, u850, v850), (500.0, u500, v500)):
        if sp_pa / 100.0 > p + 20.0 and uu == uu and vv == vv:
            z0 = R_AIR * t2m_k / G * math.log(sp_pa / (p * 100.0))
            tm = t2m_k - 0.0065 * z0 / 2.0
            nodes.append((R_AIR * tm / G * math.log(sp_pa / (p * 100.0)), uu, vv))
    if any(x != x for n in nodes[:2] for x in n):
        return float("nan"), float("nan")
    nodes.sort()
    zs = np.array([n[0] for n in nodes])
    top = max(float(blh), 10.0)
    zz = np.linspace(0.0, top, 200)
    um = float(np.mean(np.interp(zz, zs, [n[1] for n in nodes])))
    vm = float(np.mean(np.interp(zz, zs, [n[2] for n in nodes])))
    return um, vm


def _bilinear_cells(cells, slat, slon, key):
    """cells: {(lat_c, lon_c): value-dict}. Bilinear value of `key` at the site from the surrounding cell centres."""
    lats = sorted({c[0] for c in cells})
    lons = sorted({c[1] for c in cells})
    lo_la = max([x for x in lats if x <= slat + 1e-9], default=None)
    hi_la = min([x for x in lats if x >= slat - 1e-9], default=None)
    lo_lo = max([x for x in lons if x <= slon + 1e-9], default=None)
    hi_lo = min([x for x in lons if x >= slon - 1e-9], default=None)
    if None in (lo_la, hi_la, lo_lo, hi_lo):
        return float("nan")
    ty = 0.0 if hi_la == lo_la else (slat - lo_la) / (hi_la - lo_la)
    tx = 0.0 if hi_lo == lo_lo else (slon - lo_lo) / (hi_lo - lo_lo)
    g = lambda la, lo: cells[(la, lo)][key]
    return ((1 - ty) * ((1 - tx) * g(lo_la, lo_lo) + tx * g(lo_la, hi_lo))
            + ty * ((1 - tx) * g(hi_la, lo_lo) + tx * g(hi_la, hi_lo)))


def read_era5_gee(sl_csv, pl_csv, sites):
    """ERA5 exported by ERA5_GEE_export.js (CSV: 3 x 3 sample points per site and hour, 12 and 13 UTC).
    Returns the same {(site, date): record} as read_era5_set."""
    LEV = [1000, 975, 950, 925, 900, 850, 500]

    def load(path):
        grouped = {}                                   # (site, 'YYYY-MM-DD HH') -> {(lat_c, lon_c): row}
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                try:
                    k = (r["site"], r["time_utc"])
                    cell = (round(float(r["latitude"]), 4), round(float(r["longitude"]), 4))
                except (KeyError, ValueError):
                    continue
                row = {}
                for kk, vv in r.items():
                    if kk in ("site", "time_utc", "slat", "slon", "latitude", "longitude", ".geo", "system:index"):
                        continue
                    try:
                        row[kk] = float(vv)
                    except (TypeError, ValueError):
                        row[kk] = float("nan")
                grouped.setdefault(k, {})[cell] = row
        return grouped

    sl = load(sl_csv)
    pl = load(pl_csv) if pl_csv and os.path.exists(pl_csv) else {}
    out = {}
    nan = float("nan")
    for site in sites:
        la, lo = SITES[site]
        days = sorted({k[1][:10] for k in sl if k[0] == site})
        for d in days:
            hrs = [h for h in ("12", "13") if (site, "%s %s" % (d, h)) in sl]
            if not hrs:
                continue

            def at(table, h, key):
                return _bilinear_cells(table[(site, "%s %s" % (d, h))], la, lo, key)

            def mean_over(table, key):
                v = [at(table, h, key) for h in hrs if (site, "%s %s" % (d, h)) in table]
                v = [x for x in v if x == x]
                return float(np.mean(v)) if v else nan

            u10 = mean_over(sl, "u_component_of_wind_10m")
            v10 = mean_over(sl, "v_component_of_wind_10m")
            blh = mean_over(sl, "boundary_layer_height")
            t2m = mean_over(sl, "temperature_2m")
            sp = mean_over(sl, "surface_pressure")
            if any(x != x for x in (u10, v10, blh, t2m, sp)):
                continue
            rec = {"u10": u10, "v10": v10, "wspd10": math.hypot(u10, v10),
                   "wdir10": (math.degrees(math.atan2(-u10, -v10)) + 360.0) % 360.0,
                   "blh": blh, "t2m": t2m - 273.15, "sp_hpa": sp / 100.0}
            # hourly solar radiation (J/m2) in the hour ending at 13 UTC
            rec["ssrd_jm2"] = at(sl, "13", "surface_solar_radiation_downwards") if "13" in hrs else nan
            rec["tp_mm"] = nan
            rec["vc10"] = blh * rec["wspd10"]
            for k in ("wspd_layer", "wdir_layer", "n_levels_in_blh", "vclayer", "vc100", "vcprof", "wspd100", "wspd_prof", "u100", "v100", "u_prof", "v_prof"):
                rec[k] = nan
            if "u_component_of_wind_100m" in next(iter(next(iter(sl.values())).values())):
                u100 = mean_over(sl, "u_component_of_wind_100m")
                v100 = mean_over(sl, "v_component_of_wind_100m")
                if u100 == u100 and v100 == v100:
                    rec["wspd100"] = math.hypot(u100, v100)
                    rec["u100"], rec["v100"] = u100, v100
                    rec["vc100"] = blh * rec["wspd100"]
                    um, vm = profile_layer_wind(
                        u10, v10, u100, v100,
                        mean_over(sl, "u_component_of_wind_850hPa"), mean_over(sl, "v_component_of_wind_850hPa"),
                        mean_over(sl, "u_component_of_wind_500hPa"), mean_over(sl, "v_component_of_wind_500hPa"),
                        sp, t2m, blh)
                    if um == um:
                        rec["wspd_prof"] = math.hypot(um, vm)
                        rec["u_prof"], rec["v_prof"] = um, vm
                        rec["vcprof"] = blh * rec["wspd_prof"]
            if any((site, "%s %s" % (d, h)) in pl for h in hrs):
                ptop = sp / 100.0 * math.exp(-G * blh / (R_AIR * t2m))
                uu, vv = [u10], [v10]
                for p in LEV:
                    if ptop <= p <= sp / 100.0:
                        ul = mean_over(pl, "u_component_of_wind_%dhPa" % p)
                        vl = mean_over(pl, "v_component_of_wind_%dhPa" % p)
                        if ul == ul and vl == vl:
                            uu.append(ul)
                            vv.append(vl)
                ml, mv = float(np.mean(uu)), float(np.mean(vv))
                rec["wspd_layer"] = math.hypot(ml, mv)
                rec["wdir_layer"] = (math.degrees(math.atan2(-ml, -mv)) + 360.0) % 360.0
                rec["n_levels_in_blh"] = len(uu) - 1
                rec["vclayer"] = rec["blh"] * rec["wspd_layer"]
            out[(site, d)] = rec
    return out


def era5_files(era5_dir, kind, year):
    """All ERA5 files of one year and kind, whatever the cut: _YYYY, _YYYY_MM or _YYYY_mAA-BB."""
    import glob
    pat = os.path.join(era5_dir, "era5_%s_all_%d*.nc" % (kind, year))
    ok = re.compile(r"era5_%s_all_%d(_\d\d|_m\d\d-\d\d)?\.nc$" % (kind, year))
    return sorted(p for p in glob.glob(pat) if ok.search(os.path.basename(p)) and os.path.getsize(p) > 0)


def stage_merge(no2_rows, era5_dir, sites, years, out):
    met = {}
    for y in years:
        gsl = os.path.join(era5_dir, "ERA5X_SL_%d.csv" % y)     # exported from Earth Engine (with higher winds)
        if not os.path.exists(gsl):
            gsl = os.path.join(era5_dir, "ERA5_SL_%d.csv" % y)  # first Earth Engine export
        if os.path.exists(gsl):
            met.update(read_era5_gee(gsl, os.path.join(era5_dir, "ERA5_PL_%d.csv" % y), sites))
            continue
        sl = era5_files(era5_dir, "sl", y)                       # downloaded from the Copernicus CDS
        if not sl:
            continue
        met.update(read_era5_set(sl, era5_files(era5_dir, "pl", y), sites))
    print("ERA5 days read:", len(met))
    merged = []
    for r in no2_rows:
        k = (r["site"], r["date"])
        if k in met:
            rr = dict(r)
            rr.update(met[k])
            merged.append(rr)
    path = os.path.join(out, "no2_ventilation_daily.csv")
    write_csv(path, merged)
    print("wrote", path, "(%d rows with both NO2 and ERA5)" % len(merged))
    return merged


# ================================ stage 3 ====================================
def newey_west(X, resid, lag):
    n, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    Xe = X * resid[:, None]
    S = Xe.T @ Xe
    for l in range(1, lag + 1):
        w = 1.0 - l / (lag + 1.0)
        G_l = Xe[l:].T @ Xe[:-l]
        S += w * (G_l + G_l.T)
    return XtX_inv @ S @ XtX_inv


def design(dates, x_main, extra, with_month=True, with_dow=True, with_trend=True):
    dd = [dt.date.fromisoformat(d) for d in dates]
    cols = [np.ones(len(dd)), x_main] + extra
    names = ["const", "ln_vc"] + ["x%d" % i for i in range(len(extra))]
    if with_trend:
        cols.append(np.array([(d - dt.date(2019, 1, 1)).days / 365.25 for d in dd]))
        names.append("trend_per_year")
    if with_month:
        for m in range(2, 13):
            cols.append(np.array([1.0 if d.month == m else 0.0 for d in dd]))
            names.append("month%02d" % m)
    if with_dow:
        for k in range(1, 7):
            cols.append(np.array([1.0 if d.weekday() == k else 0.0 for d in dd]))
            names.append("dow%d" % k)
    return np.column_stack(cols), names


def fit_site(rows, vc_key, y_key="no2_inner", extra_keys=("t2m",), with_month=True):
    rr = [r for r in rows if r["usable"] == 1 and all(np.isfinite([r[y_key], r[vc_key]] + [r[k] for k in extra_keys]))
          and r[vc_key] > 0 and r[y_key] > 0]
    if len(rr) < 60:
        return None
    rr.sort(key=lambda r: r["date"])
    y = np.log([r[y_key] for r in rr])
    x = np.log([r[vc_key] for r in rr])
    extra = [np.array([r[k] for r in rr]) for k in extra_keys]
    X, names = design([r["date"] for r in rr], x, extra, with_month=with_month)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    V = newey_west(X, res, NW_LAG)
    se = math.sqrt(V[1, 1])
    # partial R2 of ln_vc: variance explained beyond the other regressors
    X0 = np.delete(X, 1, axis=1)
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    r0 = y - X0 @ b0
    pr2 = 1.0 - (res @ res) / (r0 @ r0)
    return {"n": len(rr), "slope_lnNO2_per_lnVC": beta[1], "se_NW": se,
            "ci95_low": beta[1] - 1.96 * se, "ci95_high": beta[1] + 1.96 * se,
            "partial_R2": pr2, "effect_of_halving_VC_pct": 100.0 * (math.exp(-beta[1] * math.log(2)) - 1.0)}


def stage_model(merged, sites, out, vc_keys=("vc10", "vc100", "vcprof", "vclayer")):
    results = []
    vc_keys = [k for k in vc_keys if sum(1 for r in merged if k in r and np.isfinite(r[k])) >= 100]
    for site in sites:
        rows = [r for r in merged if r["site"] == site]
        for vk in vc_keys:
            if not rows or vk not in rows[0]:
                continue
            fit = fit_site(rows, vk)
            if fit is None:
                print("  %s %s: too few usable days" % (site, vk))
                continue
            fit = dict(site=site, ventilation=vk, **fit)
            results.append(fit)
            print("  %-15s %-8s n=%4d  slope=%6.3f  (95%% CI %6.3f to %6.3f)  partial R2=%.3f"
                  % (site, vk, fit["n"], fit["slope_lnNO2_per_lnVC"], fit["ci95_low"], fit["ci95_high"], fit["partial_R2"]))
    write_csv(os.path.join(out, "regression_results.csv"), results)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing: no figure")
        return results
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=False)
    for ax, site in zip(axes.ravel(), sites):
        rows = [r for r in merged if r["site"] == site and r["usable"] == 1 and np.isfinite(r.get("vc10", np.nan))
                and r["no2_inner"] > 0 and r["vc10"] > 0]
        if not rows:
            ax.set_title(site + " (no data)")
            continue
        x = np.log10([r["vc10"] for r in rows])
        y = np.array([r["no2_inner"] for r in rows])
        ax.scatter(x, y, s=4, alpha=0.25, color="#888")
        q = np.quantile(x, np.linspace(0, 1, 11))
        for a, b in zip(q[:-1], q[1:]):
            m = (x >= a) & (x <= b)
            if m.sum() >= 5:
                ax.plot([x[m].mean()], [np.median(y[m])], "o", color="#c0392b")
        ax.set_yscale("log")
        ax.set_title(site)
        ax.set_xlabel("log10 ventilation (m2/s, BLH x 10 m wind)")
        ax.set_ylabel("NO2 within 10 km (umol/m2)")
    fig.suptitle("Daily NO2 against ventilation (grey: days; red: decile medians)")
    fig.tight_layout()
    fp = os.path.join(out, "no2_vs_ventilation.png")
    fig.savefig(fp, dpi=150)
    print("wrote", fp)
    return results


# ================================== main =====================================
def parse_list(text, allowed, conv=str):
    if text in (None, "", "default", "all"):
        return list(allowed)
    return [conv(t) for t in text.split(",")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=["tables", "merge", "model", "all"], default="all")
    ap.add_argument("--tif-dir", default=".")
    ap.add_argument("--era5-dir", default="era5")
    ap.add_argument("--out", default="out")
    ap.add_argument("--sites", default="all")
    ap.add_argument("--years", default="all", help="e.g. 2022 or 2019,2020")
    a = ap.parse_args()
    sites = parse_list(a.sites, DEFAULT_SITES)
    bad = [s for s in sites if s not in SITES]
    if bad:
        raise SystemExit("Unknown site(s): %s" % ", ".join(bad))
    years = [int(y) for y in parse_list(a.years, YEARS)]
    os.makedirs(a.out, exist_ok=True)

    tab = os.path.join(a.out, "no2_daily_table.csv")
    mer = os.path.join(a.out, "no2_ventilation_daily.csv")
    rows = merged = None
    if a.stage in ("tables", "all"):
        print("Stage 1: NO2 tables from GeoTIFFs")
        rows = stage_tables(a.tif_dir, sites, years, a.out)
    if a.stage in ("merge", "all"):
        if rows is None:
            rows = read_csv(tab)
        print("Stage 2: adding ERA5")
        merged = stage_merge(rows, a.era5_dir, sites, years, a.out)
    if a.stage in ("model", "all"):
        if merged is None:
            merged = read_csv(mer)
        print("Stage 3: regression")
        stage_model(merged, sites, a.out)


if __name__ == "__main__":
    main()
