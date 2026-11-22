#!/usr/bin/env python3
"""Wind-rotation analysis of TROPOMI NO2 (daily 0.01 deg GeoTIFFs from Earth Engine) with ERA5 winds.

For each usable day the pixels are rotated so that the mean boundary-layer wind blows towards +x, the rotated
fields are averaged by wind-speed class, summed across the wind (|y| <= YMAX km) into a line density, and an
exponentially modified Gaussian (Beirle et al. 2011) is fitted:  a = NO2 burden, x0 = e-folding length,
lifetime tau = x0 / w,  emission E = a / tau.  NOx emission = 1.32 * NO2 emission.

Stages:  python wind_rotation.py --stage accumulate --tif-dir DIR --era5-dir DIR --out DIR --sites algiers --years 2020
         python wind_rotation.py --stage fit --out DIR --sites algiers
"""
import argparse, os, sys, math, csv, glob, datetime as dt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_no2_ventilation as A

R_KM = 110.0                     # use pixels within this distance of the site
XMIN, XMAX, DX = -60.0, 110.0, 2.0
YMAX, DY = 40.0, 2.0
NX = int(round((XMAX - XMIN) / DX)); NY = int(round(2 * YMAX / DY))
WMIN = 2.0                        # m/s, days with weaker mean wind are dropped
CLASSES = [(2, 4), (4, 6), (6, 99), (2, 99)]
MIN_COVER = 0.5                   # share of pixels within R_KM that must be valid on a day
M_NO2 = 46.0055e-3                # kg/mol
NOX_FACTOR = 1.32
FIT_X = (-40.0, 100.0)


def read_tif(path):
    try:
        return A.read_geotiff(path)
    except SystemExit as exc:
        raise RuntimeError("GeoTIFF decoding failed; install the pinned tifffile and imagecodecs requirements.") from exc


def accumulate(tif, site, year, wind):
    arr, lon0, lat0, dx, dy = read_tif(tif)
    H, W, B = arr.shape
    slat, slon = A.SITES[site]
    lon = lon0 + (np.arange(W) + 0.5) * dx
    lat = lat0 - (np.arange(H) + 0.5) * dy
    LON, LAT = np.meshgrid(lon, lat)
    xe = (LON - slon) * 111.320 * math.cos(math.radians(slat))
    yn = (LAT - slat) * 110.574
    sel = np.hypot(xe, yn) <= R_KM
    xe, yn = xe[sel], yn[sel]
    data = arr[sel]                              # (npix, B)
    days = A.year_days(year)
    out = []
    for k, d in enumerate(days):
        rec = wind.get((site, d.isoformat()))
        if rec is None or rec.get("u_prof", float("nan")) != rec.get("u_prof", float("nan")):
            continue
        u, v = rec["u_prof"], rec["v_prof"]
        w = math.hypot(u, v)
        if w < WMIN:
            continue
        z = data[:, k]
        ok = np.isfinite(z)
        if ok.mean() < MIN_COVER:
            continue
        th = math.atan2(v, u)                    # direction the wind blows towards (math angle)
        c, s = math.cos(th), math.sin(th)
        xr = xe[ok] * c + yn[ok] * s
        yr = -xe[ok] * s + yn[ok] * c
        ix = np.floor((xr - XMIN) / DX).astype(int)
        iy = np.floor((yr + YMAX) / DY).astype(int)
        m = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        flat = ix[m] * NY + iy[m]
        S = np.bincount(flat, weights=z[ok][m], minlength=NX * NY).astype(np.float32)
        N = np.bincount(flat, minlength=NX * NY).astype(np.float32)
        out.append((d.isoformat(), w, math.degrees(th), rec["blh"], S, N))
    return out


def stage_accumulate(a):
    os.makedirs(a.out, exist_ok=True)
    for site in a.sites:
        store = {}
        for y in a.years:
            tif = os.path.join(a.tif_dir, "NO2_daily_%s_%d.tif" % (site, y))
            sl = os.path.join(a.era5_dir, "ERA5X_SL_%d.csv" % y)
            if not os.path.exists(tif):
                print("missing", tif); continue
            wind = A.read_era5_gee(sl, None, [site])
            rows = accumulate(tif, site, y, wind)
            print(site, y, len(rows), "days kept")
            for r in rows:
                store[r[0]] = r
        fn = os.path.join(a.out, "rot_%s.npz" % site)
        if os.path.exists(fn):                  # merge with earlier years
            z = np.load(fn, allow_pickle=True)
            for i, d in enumerate(z["date"]):
                store.setdefault(str(d), (str(d), z["w"][i], z["theta"][i], z["blh"][i], z["S"][i], z["N"][i]))
        ds = sorted(store)
        np.savez_compressed(fn, date=np.array(ds), w=np.array([store[d][1] for d in ds]),
                            theta=np.array([store[d][2] for d in ds]), blh=np.array([store[d][3] for d in ds]),
                            S=np.array([store[d][4] for d in ds]), N=np.array([store[d][5] for d in ds]))
        print("wrote", fn, len(ds), "days")


# ------------------------------- fitting -------------------------------------
def emg(x, a, x0, mu, sig, B):
    return (a / (2 * x0)) * np.exp(mu / x0 + sig ** 2 / (2 * x0 ** 2) - x / x0) * \
        erfc(-((x - mu) / sig - sig / x0) / math.sqrt(2)) + B


def line_density(S, N, idx):
    """S, N: (days, NX*NY). Mean field over selected days -> line density (mol/m) vs x."""
    s = S[idx].sum(0).reshape(NX, NY); n = N[idx].sum(0).reshape(NX, NY)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(n > 0, s / n, np.nan)       # mean field, umol/m2
    cols_ok = np.isfinite(f).sum(1)
    # fill missing cross-wind bins with the row mean so that every x uses the full width
    rowmean = np.nanmean(f, axis=1)
    f = np.where(np.isfinite(f), f, rowmean[:, None])
    ld = f.sum(1) * 1e-6 * (DY * 1000.0)         # umol/m2 * m  -> mol/m
    ld[cols_ok < NY * 0.6] = np.nan
    return ld


def fit_one(xc, ld, w):
    ok = np.isfinite(ld) & (xc >= FIT_X[0]) & (xc <= FIT_X[1])
    x, y = xc[ok] * 1000.0, ld[ok]               # metres
    if len(x) < 20:
        return None
    p0 = [max(y.max() - np.nanmedian(y), 1e-3) * 30e3, 40e3, 0.0, 8e3, np.nanmin(y)]
    lo = [0, 5e3, -20e3, 1e3, -np.inf]; hi = [np.inf, 400e3, 20e3, 30e3, np.inf]
    try:
        p, _ = curve_fit(emg, x, y, p0=p0, bounds=(lo, hi), maxfev=20000)
    except Exception:
        return None
    a, x0, mu, sig, B = p
    tau_h = x0 / w / 3600.0
    E = a * w / x0                                # mol/s NO2
    return dict(a_mol=a, x0_km=x0 / 1e3, mu_km=mu / 1e3, sigma_km=sig / 1e3, B_mol_m=B, tau_h=tau_h,
                E_no2_mol_s=E, E_nox_kt_yr=E * M_NO2 * NOX_FACTOR * 3.15576e7 / 1e6,
                resid_rms=float(np.sqrt(np.mean((emg(x, *p) - y) ** 2))), peak=float(y.max() - B))


def stage_fit(a):
    xc = XMIN + (np.arange(NX) + 0.5) * DX
    rng = np.random.default_rng(1)
    rows = []
    for site in a.sites:
        z = np.load(os.path.join(a.out, "rot_%s.npz" % site), allow_pickle=True)
        w, S, N = z["w"], z["S"], z["N"]
        for lo_, hi_ in CLASSES:
            idx = np.where((w >= lo_) & (w < hi_))[0]
            if len(idx) < 60:
                continue
            wm = float(np.mean(w[idx]))
            ld = line_density(S, N, idx)
            f0 = fit_one(xc, ld, wm)
            if f0 is None:
                rows.append(dict(site=site, wind_class="%g-%g" % (lo_, hi_ if hi_ < 99 else 99), n_days=len(idx), w_mean=wm, status="fit failed"))
                continue
            boots = []
            for _ in range(a.boot):
                bi = rng.choice(idx, len(idx), replace=True)
                fb = fit_one(xc, line_density(S, N, bi), float(np.mean(w[bi])))
                if fb: boots.append(fb)
            r = dict(site=site, wind_class="%g-%s" % (lo_, "up" if hi_ >= 99 else hi_), n_days=len(idx), w_mean=round(wm, 2), status="ok", n_boot=len(boots))
            for k, v in f0.items():
                r[k] = v
                if boots and k in ("tau_h", "E_nox_kt_yr", "x0_km", "a_mol"):
                    q = np.percentile([b[k] for b in boots], [2.5, 97.5])
                    r[k + "_lo"], r[k + "_hi"] = q[0], q[1]
            rows.append(r)
            np.save(os.path.join(a.out, "ld_%s_%s.npy" % (site, r["wind_class"])), np.vstack([xc, ld]))
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    fn = os.path.join(a.out, "wind_rotation_results.csv")
    with open(fn, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r in rows: wr.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k in ("site", "wind_class", "n_days", "w_mean", "x0_km", "tau_h", "tau_h_lo", "tau_h_hi", "E_nox_kt_yr", "sigma_km", "mu_km", "peak", "status")})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["accumulate", "fit"], required=True)
    ap.add_argument("--tif-dir"); ap.add_argument("--era5-dir"); ap.add_argument("--out", required=True)
    ap.add_argument("--sites", nargs="+", default=list(A.DEFAULT_SITES)); ap.add_argument("--years", nargs="+", type=int, default=list(range(2019, 2026)))
    ap.add_argument("--boot", type=int, default=100)
    a = ap.parse_args()
    stage_accumulate(a) if a.stage == "accumulate" else stage_fit(a)
