#!/usr/bin/env python3
"""Figure 1: daily NO2 (10 km disc) against ventilation VC10 at the six study sites.
Usage: python fig_no2_vs_ventilation.py --table no2_ventilation_daily_v2.csv --out fig1_no2_vs_ventilation.png
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_no2_ventilation as A

SITES = [("algiers", "Algiers"), ("oran", "Oran"), ("arzew", "Arzew"),
         ("constantine", "Constantine"), ("setif", "Sétif"), ("hassi_messaoud", "Hassi Messaoud")]
INK, GREY, RED = "#2b2b2b", "#9a9a95", "#b5421f"
ap = argparse.ArgumentParser(); ap.add_argument("--table", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = A.read_csv(a.table)
fig, axes = plt.subplots(2, 3, figsize=(11, 6.6), dpi=200)
for ax, (key, name) in zip(axes.ravel(), SITES):
    rows_for_beta = [r for r in rows if r["site"] == key]
    rr = [r for r in rows if r["site"] == key and r["usable"] == 1 and r["vc10"] > 0 and r["no2_inner"] > 0
          and np.isfinite(r["vc10"])]
    x = np.array([r["vc10"] for r in rr]); y = np.array([r["no2_inner"] for r in rr])
    ax.scatter(x, y, s=5, color=GREY, alpha=0.3, linewidths=0, zorder=2)
    q = np.quantile(np.log10(x), np.linspace(0, 1, 11))
    mx, my = [], []
    for lo, hi in zip(q[:-1], q[1:]):
        m = (np.log10(x) >= lo) & (np.log10(x) <= hi)
        if m.sum() >= 5:
            mx.append(10 ** np.log10(x[m]).mean()); my.append(np.median(y[m]))
    ax.plot(mx, my, "-o", color=RED, ms=5, lw=1.6, mec="#fcfcfb", mew=0.8, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title("%s (n = %d)" % (name, len(rr)), fontsize=10.5, color=INK, loc="left")
    ax.text(0.97, 0.95, "β = %.2f" % A.fit_site(rows_for_beta, "vc10")["slope_lnNO2_per_lnVC"], transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=INK)
    ax.grid(color="#e7e7e3", lw=0.7, zorder=0); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9, colors=INK)
for i, ax in enumerate(axes.ravel()):
    if i % 3 == 0: ax.set_ylabel("NO₂ in 10 km disc (µmol m⁻²)", fontsize=9.5, color=INK)
    if i >= 3: ax.set_xlabel("Ventilation VC10 (m² s⁻¹)", fontsize=9.5, color=INK)
fig.tight_layout()
fig.savefig(a.out, facecolor="#fcfcfb")
