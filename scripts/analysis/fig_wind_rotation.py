#!/usr/bin/env python3
"""Figure S1: wind-rotated NO2 line densities and exponentially modified Gaussian fits (all days, wind >= 2 m/s).

Reads outputs/wind_rotation/ld_<site>_2-up.npy (x in km, line density in mol/m) and
outputs/wind_rotation/wind_rotation_results.csv (fit parameters), both written by scripts/run/3_wind_rotation.bat.
Writes outputs/fig_wind_rotation.png.
Usage:  python scripts\\analysis\\fig_wind_rotation.py
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--input-dir',required=True);ap.add_argument('--out',required=True);args=ap.parse_args()
DIR=args.input_dir
OUT=args.out

SITES = [("algiers", "Algiers", False), ("oran", "Oran", False), ("arzew", "Arzew", False),
         ("constantine", "Constantine", True), ("setif", "Sétif", True), ("hassi_messaoud", "Hassi Messaoud", False)]
FIT_X = (-40.0, 100.0)


def emg(x, a, x0, mu, sig, B):
    return (a / (2 * x0)) * np.exp(mu / x0 + sig ** 2 / (2 * x0 ** 2) - x / x0) * \
        erfc(-((x - mu) / sig - sig / x0) / math.sqrt(2)) + B


def main():
    rows = {}
    with open(os.path.join(DIR, "wind_rotation_results.csv")) as f:
        for r in csv.DictReader(f):
            if r["wind_class"] == "2-up" and r["status"] == "ok":
                rows[r["site"]] = r
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.6))
    for ax, (key, name, bad) in zip(axes.ravel(), SITES):
        xc, ld = np.load(os.path.join(DIR, "ld_%s_2-up.npy" % key))
        r = rows[key]
        ax.plot(xc, ld, "o", ms=3.2, color="#6b6b6b", label="Wind-rotated mean", zorder=2)
        xs = np.linspace(FIT_X[0], FIT_X[1], 300)
        ys = emg(xs * 1000.0, float(r["a_mol"]), float(r["x0_km"]) * 1000.0, float(r["mu_km"]) * 1000.0,
                 float(r["sigma_km"]) * 1000.0, float(r["B_mol_m"]))
        ax.plot(xs, ys, lw=2.0, color="#c0392b" if not bad else "#d98c84",
                ls="-" if not bad else "--", label="Exponentially modified Gaussian fit", zorder=3)
        tau, lo, hi = float(r["tau_h"]), float(r["tau_h_lo"]), float(r["tau_h_hi"])
        e = float(r["E_nox_kt_yr"])
        title = name + (" (fit not reliable)" if bad else "")
        ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left")
        ax.text(0.97, 0.95, "τ = %.1f h (%.1f–%.1f)\nE(NOₓ) = %.0f kt yr⁻¹\n%d days" % (tau, lo, hi, e, int(r["n_days"])),
                transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color="#333333",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.6))
        ax.set_xlabel("Distance downwind (km)", fontsize=9)
        ax.set_ylabel("NO₂ line density (mol m⁻¹)", fontsize=9)
        ax.tick_params(labelsize=8.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(alpha=0.25, lw=0.5)
        # leave room for the text box without hiding the plume
        y0, y1 = ax.get_ylim()
        ax.set_ylim(y0, y1 + 0.18 * (y1 - y0))
    h, l = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Wind-rotated NO₂ line densities, all days with ERA5 wind ≥ 2 m s⁻¹, 2019–2025", fontsize=11.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(OUT, dpi=170)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
