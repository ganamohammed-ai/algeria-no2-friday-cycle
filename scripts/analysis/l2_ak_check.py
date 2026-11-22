#!/usr/bin/env python3
"""Level-2 check with averaging kernels: how much does the assumed NO2 profile shape matter?

For every orbit clip (data\\02_tropomi_L2_orbit_clips\\<site>\\*.npz) and the pixels inside the
10 km disc around the site (quality value >= 0.75) it
  1. builds the tropospheric averaging kernel  A_k = AK_k * AMF_total / AMF_trop  (layers up to the TM5 tropopause),
  2. rebuilds the TM5 layer pressures from tm5_constant_a/b and the surface pressure,
  3. replaces the TM5 a-priori profile by two simple alternatives, using the ERA5 boundary-layer height of that day:
       "bl"      : all tropospheric NO2 mixed uniformly from the surface to the top of the boundary layer
       "uniform" : NO2 mixing ratio uniform through the whole troposphere
  4. recomputes the column:  VCD_new = VCD_old / ( sum_k A_k x_k / sum_k x_k ),  x_k = layer partial column of the new profile.
It writes one row per site and day (means over the pixels), small enough to send back.

Usage (from the project folder):
  python scripts\\analysis\\l2_ak_check.py --clips data\\02_tropomi_L2_orbit_clips --table outputs\\no2_ventilation_daily_v2.csv --out outputs
"""
import argparse, csv, glob, math, os, sys, collections
import numpy as np

np.seterr(divide="ignore", invalid="ignore")

SITES = {"algiers": (36.75, 3.06), "oran": (35.70, -0.63), "arzew": (35.85, -0.31),
         "constantine": (36.37, 6.61), "setif": (36.19, 5.41), "hassi_messaoud": (31.68, 6.07)}
RADIUS_KM = 10.0
QA_MIN = 0.75
R_EARTH = 6371.0
RD, G = 287.05, 9.80665


def layer_pressures(a, b, ps):
    """a, b: (34,2) [bottom, top] interface coefficients; ps: (n,) Pa. Returns (n,34) bottom and top pressure."""
    pb = a[None, :, 0] + b[None, :, 0] * ps[:, None]
    pt = a[None, :, 1] + b[None, :, 1] * ps[:, None]
    return pb, pt


def new_shape_weights(pb, pt, tidx, p_bl):
    """Partial-column weights x_k (n,34) for the two alternative profiles (zero above the tropopause layer)."""
    k = np.arange(pb.shape[1])[None, :]
    trop = k <= tidx[:, None]
    thick = np.where(trop, np.clip(pb - pt, 0, None), 0.0)                 # pressure thickness in the troposphere
    in_bl = np.clip(pb - np.maximum(pt, p_bl[:, None]), 0, None)            # overlap with [p_bl, ps]
    x_bl = np.where(trop, np.minimum(in_bl, thick), 0.0)
    return x_bl, thick


def compute_file(z, lat0, lon0, blh_m, t2m_c):
    """Return dict of pixel sums for one orbit clip, or None."""
    lat = z["latitude"].astype(float); lon = z["longitude"].astype(float)
    dx = np.radians(lon - lon0) * R_EARTH * math.cos(math.radians(lat0))
    dy = np.radians(lat - lat0) * R_EARTH
    ok = (dx ** 2 + dy ** 2 <= RADIUS_KM ** 2) & (z["qa_value"] >= QA_MIN) & np.isfinite(z["no2_trop_mol_m2"])
    ok &= np.isfinite(z["air_mass_factor_troposphere"]) & (z["air_mass_factor_troposphere"] > 0)
    if not ok.any():
        return None
    s = lambda name: z[name][ok].astype(float)
    ak = s("averaging_kernel"); mt = s("air_mass_factor_total"); mtr = s("air_mass_factor_troposphere")
    ps = s("surface_pressure"); tidx = z["tm5_tropopause_layer_index"][ok].astype(int)
    vcd = s("no2_trop_mol_m2") * 1e6                                         # umol m-2, like the main table
    atrop = ak * (mt / mtr)[:, None]
    a = z["tm5_constant_a"].astype(float); b = z["tm5_constant_b"].astype(float)
    pb, pt = layer_pressures(a, b, ps)
    H = RD * (t2m_c + 273.15 - 3.0) / G                                     # scale height, m
    p_bl = ps * math.exp(-blh_m / H)                                        # (n,) Pa
    x_bl, x_un = new_shape_weights(pb, pt, tidx, p_bl)
    def factor(x):
        den = x.sum(axis=1)
        return np.where(den > 0, (atrop * x).sum(axis=1) / np.where(den > 0, den, 1), np.nan)
    f_bl, f_un = factor(x_bl), factor(x_un)
    good = np.isfinite(f_bl) & np.isfinite(f_un) & (f_bl > 0.05) & (f_un > 0.05)
    if not good.any():
        return None
    g = good
    t = z["time_utc"][ok][g]
    hrs = np.array([int(u[11:13]) + int(u[14:16]) / 60.0 for u in t])
    return {"n": int(g.sum()),
            "vcd": vcd[g].sum(), "vcd_bl": (vcd / f_bl)[g].sum(), "vcd_un": (vcd / f_un)[g].sum(),
            "f_bl": f_bl[g].sum(), "f_un": f_un[g].sum(),
            "amf": mtr[g].sum(), "cf": s("cloud_radiance_fraction_nitrogendioxide_window")[g].sum(),
            "alb": s("surface_albedo_nitrogendioxide_window")[g].sum(), "hour": hrs.sum(),
            "bl_frac_layers": (x_bl[g] > 0).sum()}


def read_blh(path):
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                out[(r["site"], r["date"])] = (float(r["blh"]), float(r["t2m"]))
            except (ValueError, KeyError):
                pass
    return out


def main():
    global RADIUS_KM, QA_MIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--clips", required=True); ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--radius", type=float, default=RADIUS_KM); ap.add_argument("--qa", type=float, default=QA_MIN)
    a = ap.parse_args()
    RADIUS_KM, QA_MIN = a.radius, a.qa
    meteo = read_blh(a.table)
    rows = []
    skipped = collections.Counter()
    for site, (la, lo) in SITES.items():
        files = sorted(glob.glob(os.path.join(a.clips, site, "*.npz")))
        agg = {}
        for i, fpath in enumerate(files):
            d8 = os.path.basename(fpath)[:8]; date = "%s-%s-%s" % (d8[:4], d8[4:6], d8[6:8])
            m = meteo.get((site, date))
            if m is None:
                skipped["no_meteo"] += 1; continue
            try:
                z = np.load(fpath, allow_pickle=True)
                res = compute_file(z, la, lo, m[0], m[1])
            except Exception as e:                                       # a damaged clip must not stop the run
                skipped["error"] += 1; print("skip", fpath, e); continue
            if res is None:
                skipped["no_pixels"] += 1; continue
            cur = agg.setdefault(date, collections.Counter())
            for k, v in res.items():
                cur[k] += v
            if (i + 1) % 200 == 0:
                print(site, i + 1, "/", len(files), flush=True)
        for date in sorted(agg):
            c = agg[date]; n = c["n"]
            rows.append({"site": site, "date": date, "n_pix": n,
                         "vcd": c["vcd"] / n, "vcd_bl": c["vcd_bl"] / n, "vcd_un": c["vcd_un"] / n,
                         "f_bl": c["f_bl"] / n, "f_un": c["f_un"] / n, "amf_trop": c["amf"] / n,
                         "cloud_rad_frac": c["cf"] / n, "albedo": c["alb"] / n, "hour_utc": c["hour"] / n})
        print(site, "days:", len(agg))
    out = os.path.join(a.out, "l2_ak_daily.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("wrote", out, "rows", len(rows), "skipped", dict(skipped))


if __name__ == "__main__":
    main()
