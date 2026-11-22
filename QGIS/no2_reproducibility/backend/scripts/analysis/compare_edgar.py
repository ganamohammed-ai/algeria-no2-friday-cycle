"""Compare wind-rotation NOx emissions with a gridded inventory (EDGAR or CAMS-GLOB-ANT).

Run on your own computer after downloading ONE NetCDF file with annual total NOx
(EDGAR v8.1: https://edgar.jrc.ec.europa.eu/dataset_ap81 , NOx, year 2022 or the nearest year,
the annual TOTALS "flx" NetCDF, kg m-2 s-1; an "emi" file in tonnes per cell also works).
Needs: pip install netCDF4 numpy

Usage (from scripts\\analysis):
  python compare_edgar.py --nc C:\\path\\EDGAR_NOx_2022.nc --rotation ..\\..\\outputs\\wind_rotation_results.csv --out ..\\..\\outputs

For each site it sums the inventory inside discs of 20, 40 and 80 km radius and writes
inventory_comparison.csv next to the wind-rotation value (class 2-up).
The comparison assumes compatible NOx mass conventions. Verify the inventory
metadata before interpreting ratios; the manuscript does not claim that this
mass-basis check or a matching satellite footprint was established.
Discs are cut at 0.1 deg cell centres, so the 20 km sum is about 10-15% rough; use the 40 and 80 km sums for comparison.
The inventory year is NOT the same as the satellite period (2019-2025): say so in the paper.
"""
import argparse, csv, math, os
import numpy as np

SITES = {"algiers": (36.75, 3.06), "oran": (35.70, -0.63), "arzew": (35.85, -0.31),
         "constantine": (36.37, 6.61), "setif": (36.19, 5.41), "hassi_messaoud": (31.68, 6.07)}
RADII_KM = (20, 40, 80)
R_EARTH = 6371.0


def cell_area_km2(lat, dlat, dlon):
    """Area of cells centred on lat (deg, array) with size dlat x dlon (deg)."""
    return (R_EARTH ** 2) * math.radians(dlon) * (
        np.sin(np.radians(lat + dlat / 2)) - np.sin(np.radians(lat - dlat / 2)))


def disc_sum(lat, lon, grid, lat0, lon0, radius_km, per_area):
    """Sum of grid inside a disc. per_area=True: grid is a flux per m2 per s (-> multiplied by area).
    per_area=False: grid is mass per cell. Returns the same unit as the grid (per cell) summed."""
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    dlat = abs(lat[1] - lat[0]); dlon = abs(lon[1] - lon[0])
    LON, LAT = np.meshgrid(lon, lat)
    dx = np.radians(LON - lon0) * R_EARTH * math.cos(math.radians(lat0))
    dy = np.radians(LAT - lat0) * R_EARTH
    inside = (dx ** 2 + dy ** 2) <= radius_km ** 2
    g = np.where(np.isfinite(grid), grid, 0.0)
    if per_area:
        g = g * cell_area_km2(LAT, dlat, dlon) * 1e6
    return float(g[inside].sum()), int(inside.sum())


def read_nc(path):
    from netCDF4 import Dataset
    ds = Dataset(path)
    names = list(ds.variables)
    latn = next(n for n in names if n.lower() in ("lat", "latitude"))
    lonn = next(n for n in names if n.lower() in ("lon", "longitude"))
    cand = [n for n in names if n not in (latn, lonn) and ds.variables[n].ndim >= 2
            and not n.lower().startswith(("time", "crs", "bnds"))]
    if not cand:
        raise SystemExit("no 2-D variable found in " + path)
    v = ds.variables[cand[0]]
    units = getattr(v, "units", "")
    arr = np.ma.filled(v[:], np.nan).astype(float)
    while arr.ndim > 2:
        arr = arr[0]
    lat = np.array(ds.variables[latn][:]); lon = np.array(ds.variables[lonn][:])
    if lon.max() > 180:
        lon = np.where(lon > 180, lon - 360, lon)
        order = np.argsort(lon); lon = lon[order]; arr = arr[:, order]
    if lat[0] > lat[-1]:
        lat = lat[::-1]; arr = arr[::-1, :]
    return cand[0], units, lat, lon, arr


def to_kt_per_year(total, units):
    """Convert a disc sum to kt/yr from the grid's unit string."""
    u = units.replace(" ", "").lower()
    if "kgm-2s-1" in u or "kg/m2/s" in u or "kgm^-2s^-1" in u:
        return total * 365.25 * 86400 / 1e6          # total is kg/s after area multiplication
    if "tonnes" in u or "ton" in u or u in ("t", "t/year", "tonne/year"):
        return total / 1e3                              # tonnes/yr -> kt/yr
    if u in ("kg", "kg/year", "kgyr-1"):
        return total / 1e6
    raise SystemExit("unknown unit %r: convert by hand" % units)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nc", required=True)
    ap.add_argument("--rotation", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    var, units, lat, lon, arr = read_nc(a.nc)
    per_area = any(k in units.replace(" ", "").lower() for k in ("m-2", "/m2"))
    print("variable %s, units %r, per-area=%s" % (var, units, per_area))
    rot = {}
    for r in csv.DictReader(open(a.rotation)):
        if r["wind_class"] == "2-up":
            rot[r["site"]] = r
    rows = []
    for s, (la, lo) in SITES.items():
        row = {"site": s}
        for R in RADII_KM:
            tot, ncell = disc_sum(lat, lon, arr, la, lo, R, per_area)
            row["inventory_%dkm_kt_yr" % R] = round(to_kt_per_year(tot, units), 2)
        r = rot.get(s)
        row["rotation_kt_yr"] = r["E_nox_kt_yr"] if r else ""
        row["rotation_lo"] = r["E_nox_kt_yr_lo"] if r else ""
        row["rotation_hi"] = r["E_nox_kt_yr_hi"] if r else ""
        rows.append(row)
        print(row)
    out = os.path.join(a.out, "inventory_comparison.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
