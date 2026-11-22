#!/usr/bin/env python3
"""Figure 4: change in NO2 during 23 Mar-13 Jun 2020, raw vs weather-adjusted, with 95% CIs.
Usage: python fig_lockdown.py --effects lockdown_effect.csv --out fig_lockdown_effect.png"""
import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAMES = {"algiers": "Algiers", "oran": "Oran", "arzew": "Arzew", "constantine": "Constantine",
         "setif": "Sétif", "hassi_messaoud": "Hassi Messaoud",
         "control_tell": "Control (Tell)", "control_highland": "Control (highland)"}
ap = argparse.ArgumentParser()
ap.add_argument("--effects", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = [r for r in csv.DictReader(open(a.effects)) if r["window"] == "primary_0323_0613"]
order = list(NAMES)
rows.sort(key=lambda r: order.index(r["site"]))
BLUE, ORANGE, INK, GRID = "#1f5fa8", "#c2591a", "#2b2b2b", "#e3e3df"
fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=200)
ypos = []
y = 0
for i, r in enumerate(rows):
    if i == 6:
        y += 0.6
    ypos.append(y)
    y += 1
for r, y in zip(rows, ypos):
    for key, col, mk, dy, lab in (("raw", BLUE, "o", -0.14, "Unadjusted"), ("adj", ORANGE, "s", 0.14, "Weather-adjusted")):
        v, lo, hi = float(r[key + "_pct"]), float(r[key + "_lo"]), float(r[key + "_hi"])
        ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=2, solid_capstyle="round", zorder=2)
        ax.plot([v], [y + dy], marker=mk, color=col, ms=7, mec="#fcfcfb", mew=1.2, zorder=3,
                label=lab if r is rows[0] else None)
ax.axvline(0, color=INK, lw=1, zorder=1)
ax.axhline(ypos[5] + 0.8, color=GRID, lw=1)
ax.set_yticks(ypos)
ax.set_yticklabels([NAMES[r["site"]] for r in rows], fontsize=10, color=INK)
ax.invert_yaxis()
ax.set_xlabel("Change in NO₂ column, 23 March to 13 June 2020 (%, 95% CI)", fontsize=10, color=INK)
ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
ax.set_axisbelow(True)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=10)
fig.tight_layout()
fig.savefig(a.out, facecolor="#fcfcfb")
